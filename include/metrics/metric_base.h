#ifndef _METRIC_BASE_H_
#define _METRIC_BASE_H_

#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <immintrin.h>
#include <stdexcept>
#include "CudaLE.h"
#include "vec3.h"
#include "constant_defs.h"
#include "utils/util_functions.h"
#include "metrics.h"

namespace Aperture {

namespace metric {

using CudaLE::placeholders::_1;
using CudaLE::placeholders::_2;
using CudaLE::placeholders::_3;
using CudaLE::pow;
using CudaLE::ZeroOp;

struct metric_base_interface {
  virtual double operator()(int i, int j, double x1, double x2,
                            double x3) const = 0;
  virtual double operator()(int i, int j, double x[]) const = 0;
  virtual double alpha(double x1, double x2, double x3) const = 0;
  virtual double alpha(double x[]) const = 0;
  virtual double beta(int i, double x1, double x2, double x3) const = 0;
  virtual double beta(int i, double x[]) const = 0;

  virtual double det(double x1, double x2, double x3) const = 0;
  virtual double det(double x[]) const = 0;

  virtual MetricType type() const = 0;
  virtual bool is_GR() const = 0;
  virtual bool is_linear() const = 0;
  virtual bool swap_zy() const = 0;
  virtual bool is_diagonal() const = 0;
};

template <class Impl>
struct metric_base : public metric_base_interface {
  enum { isLinear = 0 };
  enum { swapZY = 0 };
  enum { isGR = 0 };
  enum { isDiagonal = 1 };
  enum { canMapCartesian = 1 };

  double operator()(int i, int j, double x1, double x2, double x3) const {
    if (i == 1 && j == 1) {
      return static_cast<const Impl*>(this)->g11(x1, x2, x3);
    } else if ((i == 1 && j == 2) || (i == 2 && j == 1)) {
      return static_cast<const Impl*>(this)->g12(x1, x2, x3);
    } else if ((i == 1 && j == 3) || (i == 3 && j == 1)) {
      return static_cast<const Impl*>(this)->g13(x1, x2, x3);
    } else if (i == 2 && j == 2) {
      return static_cast<const Impl*>(this)->g22(x1, x2, x3);
    } else if ((i == 2 && j == 3) || (i == 3 && j == 2)) {
      return static_cast<const Impl*>(this)->g23(x1, x2, x3);
    } else if (i == 3 && j == 3) {
      return static_cast<const Impl*>(this)->g33(x1, x2, x3);
    } else if (i == 0 && j == 0) {
      return static_cast<const Impl*>(this)->g00(x1, x2, x3);
    } else if ((i == 0 && j == 1) || (i == 1 && j == 0)) {
      return static_cast<const Impl*>(this)->g01(x1, x2, x3);
    } else if ((i == 0 && j == 2) || (i == 2 && j == 0)) {
      return static_cast<const Impl*>(this)->g02(x1, x2, x3);
    } else if ((i == 0 && j == 3) || (i == 3 && j == 0)) {
      return static_cast<const Impl*>(this)->g02(x1, x2, x3);
    } else {
      throw std::runtime_error("Unknown component of the metric!");
    }
  }

  double operator()(int i, int j, double x[]) const {
    return operator()(i, j, x[0], x[1], x[2]);
  }

  double alpha(double x1, double x2, double x3) const {
    return static_cast<const Impl*>(this)->a(x1, x2, x3);
  }

  double alpha(double x[]) const { return alpha(x[0], x[1], x[2]); }

  double beta(int i, double x1, double x2, double x3) const {
    if (i == 1) {
      return static_cast<const Impl*>(this)->b1(x1, x2, x3);
    } else if (i == 2) {
      return static_cast<const Impl*>(this)->b2(x1, x2, x3);
    } else if (i == 3) {
      return static_cast<const Impl*>(this)->b3(x1, x2, x3);
    } else {
      return 0.0;
    }
  }

  double beta(int i, double x[]) const { return beta(i, x[0], x[1], x[2]); }

  double det(double x1, double x2, double x3) const {
    auto g11 = static_cast<const Impl*>(this)->g11;
    auto g22 = static_cast<const Impl*>(this)->g22;
    auto g33 = static_cast<const Impl*>(this)->g33;
    auto g12 = static_cast<const Impl*>(this)->g12;
    auto g13 = static_cast<const Impl*>(this)->g13;
    auto g23 = static_cast<const Impl*>(this)->g23;
    return (g11 * g22 * g33 - g11 * g23 * g23 + g12 * g23 * g13 -
                g12 * g12 * g33 + g13 * g12 * g23 -
                g13 * g22 * g13)(x1, x2, x3);
  }

  double det(double x[]) const { return det(x[0], x[1], x[2]); }

  MetricType type() const { return (MetricType)Impl::type; }
  bool is_GR() const { return Impl::isGR; }
  bool is_linear() const { return Impl::isLinear; }
  bool swap_zy() const { return Impl::swapZY; }
  bool is_diagonal() const { return Impl::isDiagonal; }

  DEFINE_FUNCTOR(g00, CudaLE::ConstOp(-1.0));

  DEFINE_FUNCTOR(g01, ZeroOp{});
  DEFINE_FUNCTOR(g02, ZeroOp{});
  DEFINE_FUNCTOR(g03, ZeroOp{});
  DEFINE_FUNCTOR(g10, ZeroOp{});
  DEFINE_FUNCTOR(g20, ZeroOp{});
  DEFINE_FUNCTOR(g30, ZeroOp{});

  DEFINE_FUNCTOR(g11, CudaLE::ConstOp(1.0));
  DEFINE_FUNCTOR(g22, CudaLE::ConstOp(1.0));
  DEFINE_FUNCTOR(g33, CudaLE::ConstOp(1.0));

  DEFINE_FUNCTOR(g12, ZeroOp{});
  DEFINE_FUNCTOR(g21, ZeroOp{});
  DEFINE_FUNCTOR(g13, ZeroOp{});
  DEFINE_FUNCTOR(g31, ZeroOp{});
  DEFINE_FUNCTOR(g23, ZeroOp{});
  DEFINE_FUNCTOR(g32, ZeroOp{});

  DEFINE_FUNCTOR(a, CudaLE::ConstOp(1.0));
  DEFINE_FUNCTOR(b1, ZeroOp{});
  DEFINE_FUNCTOR(b2, ZeroOp{});
  DEFINE_FUNCTOR(b3, ZeroOp{});
};

}

namespace detail {

// TODO: optimize these functions a bit!!

// template <typename F, typename G>
// general_sph_pos_to_cart(Vec3<Scalar>& pos, const F& radius, const G& g33) {
// HD_INLINE void
HOST_DEVICE void
general_sph_pos_to_cart(Vec3<Scalar>& pos);

template <typename F, typename G>
HD_INLINE void
general_sph_pos_to_cart(Vec3<Scalar>& pos, const F& radius, const G& g33) {
  double x1n = pos.x, x2n = pos.y, x3n = pos.z;
  // These become x, y, z
  pos.x = std::sqrt(g33(x1n, x2n)) * cos(x3n);
  pos.y = std::sqrt(g33(x1n, x2n)) * sin(x3n);
  // pos.z = sgn(CONST_PI * 0.5 - x2n) * std::sqrt(g11(x1n, x2n) - g33(x1n, x2n));
  pos.z = radius(x1n, x2n) * cos(x2n);
}

HOST_DEVICE void
general_sph_pos_from_cart(Vec3<Scalar>& pos);

HOST_DEVICE void
general_sph_vec_to_cart(Vec3<Scalar>& v, const Vec3<Scalar>& pos);

HOST_DEVICE void
general_sph_vec_from_cart(Vec3<Scalar>& v, const Vec3<Scalar>& pos);

HOST_DEVICE void
log_sph_pos_to_cart(Vec3<Scalar>& pos);

HOST_DEVICE void
log_sph_pos_from_cart(Vec3<Scalar>& pos);

#if defined(__AVX2__) && (defined(__ICC) || defined(__INTEL_COMPILER))
void
log_sph_pos_to_cart(__m256d* pos_x, __m256d* pos_y, __m256d* pos_z);

void
log_sph_pos_from_cart(__m256d* pos_x, __m256d* pos_y, __m256d* pos_z);

void
general_sph_vec_to_cart(__m256d* v1, __m256d* v2, __m256d* v3,
                        const __m256d* pos_x, const __m256d* pos_y, const __m256d* pos_z);

void
general_sph_vec_from_cart(__m256d* v1, __m256d* v2, __m256d* v3,
                          const __m256d* pos_x, const __m256d* pos_y, const __m256d* pos_z);
#endif

}


}

#endif  // _METRIC_BASE_H_
