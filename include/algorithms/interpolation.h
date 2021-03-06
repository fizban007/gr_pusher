#ifndef _INTERPOLATION_H_
#define _INTERPOLATION_H_

#include "cuda/cuda_control.h"
#include "typedefs.h"
#include "vec3.h"
#include <cmath>
#include <utility>
// #include "boost/fusion/container/vector.hpp"

#ifndef __CUDACC__
using std::abs;
using std::max;
#endif

namespace Aperture {

namespace detail {

struct interp_nearest_grid_point {
  enum { radius = 1, support = 1 };

  HD_INLINE Scalar operator()(Scalar dx) const {
    return (std::abs(dx) <= 0.5 ? 1.0 : 0.0);
  }
};

struct interp_cloud_in_cell {
  enum { radius = 1, support = 2 };

  HD_INLINE Scalar operator()(Scalar dx) const {
    return max(1.0 - std::abs(dx), 0.0);
  }
};

struct interp_triangular_shaped_cloud {
  enum { radius = 2, support = 3 };

  HD_INLINE Scalar operator()(Scalar dx) const {
    Scalar abs_dx = std::abs(dx);
    if (abs_dx < 0.5) {
      return 0.75 - dx * dx;
    } else if (abs_dx < 1.5) {
      Scalar tmp = 1.5 - abs_dx;
      return 0.5 * tmp * tmp;
    } else {
      return 0.0;
    }
  }
};

struct interp_piecewise_cubic {
  enum { radius = 2, support = 4 };

  HD_INLINE Scalar operator()(Scalar dx) const {
    Scalar abs_dx = std::abs(dx);
    if (abs_dx < 1.0) {
      Scalar tmp = abs_dx * abs_dx;
      return 2.0 / 3.0 - tmp + 0.5 * tmp * abs_dx;
    } else if (abs_dx < 2.0) {
      Scalar tmp = 2.0 - abs_dx;
      return 1.0 / 6.0 * tmp * tmp * tmp;
    } else {
      return 0.0;
    }
  }
};
}

// template <int Order>
class Interpolator {
 public:
  HOST_DEVICE Interpolator(int order = 1) : m_order(order) {}
  HOST_DEVICE ~Interpolator() {}

  template <typename FloatT>
  HD_INLINE Scalar interp_cell(FloatT pos, int p_cell, int target_cell,
                               int stagger = 0) const {
    FloatT x =
        (Scalar)target_cell + (stagger == 0 ? 0.5 : 1.0) - pos - (Scalar)p_cell;
    switch (m_order) {
      case 0:
        return m_interp_0(x);
      case 1:
        return m_interp_1(x);
      case 2:
        return m_interp_2(x);
      case 3:
        return m_interp_3(x);
      default:
        return 0.0;
    }
  }

  template <typename FloatT>
  HD_INLINE Vec3<Scalar> interp_cell(
      const Vec3<FloatT>& pos, const Vec3<int>& p_cell,
      const Vec3<int>& target_cell,
      Stagger_t stagger = Stagger_t("000")) const {
    return Vec3<Scalar>(
        interp_cell(pos[0], p_cell[0], target_cell[0], stagger[0]),
        interp_cell(pos[1], p_cell[1], target_cell[1], stagger[1]),
        interp_cell(pos[2], p_cell[2], target_cell[2], stagger[2]));
  }

  template <typename FloatT>
  HD_INLINE Scalar interp_weight(const Vec3<FloatT>& pos,
                                 const Vec3<int>& p_cell,
                                 const Vec3<int>& target_cell,
                                 Stagger_t stagger = Stagger_t("000")) const {
    return interp_cell(pos[0], p_cell[0], target_cell[0], stagger[0]) *
           interp_cell(pos[1], p_cell[1], target_cell[1], stagger[1]) *
           interp_cell(pos[2], p_cell[2], target_cell[2], stagger[2]);
  }

  HD_INLINE int radius() const {
    switch (m_order) {
      case 0:
        return detail::interp_nearest_grid_point::radius;
      case 1:
        return detail::interp_cloud_in_cell::radius;
      case 2:
        return detail::interp_triangular_shaped_cloud::radius;
      case 3:
        return detail::interp_piecewise_cubic::radius;
      default:
        return 0.0;
    }
  }

  HD_INLINE int support() const {
    switch (m_order) {
      case 0:
        return detail::interp_nearest_grid_point::support;
      case 1:
        return detail::interp_cloud_in_cell::support;
      case 2:
        return detail::interp_triangular_shaped_cloud::support;
      case 3:
        return detail::interp_piecewise_cubic::support;
      default:
        return 0.0;
    }
  }

 private:
  detail::interp_nearest_grid_point m_interp_0;
  detail::interp_cloud_in_cell m_interp_1;
  detail::interp_triangular_shaped_cloud m_interp_2;
  detail::interp_piecewise_cubic m_interp_3;
  int m_order;
};  // ----- end of class interpolator -----
}

#endif  // _INTERPOLATION_H_
