#ifndef _METRIC_LOG_SPHERICAL_H_
#define _METRIC_LOG_SPHERICAL_H_

#include "metric_base.h"

namespace Aperture {

namespace metric {

using CudaLE::placeholders::spherical::_r;
using CudaLE::placeholders::spherical::_theta;
using CudaLE::placeholders::spherical::_phi;

struct metric_log_spherical : public metric_base<metric_log_spherical> {
  enum { isLinear = 0 };
  enum { swapZY = 1 };
  enum { isGR = 0 };
  enum { isDiagonal = 1 };
  enum { canMapCartesian = 1 };
  enum { type = (int)MetricType::log_spherical };

  std::string name() const { return "Log Spherical"; }

  DEFINE_FUNCTOR(g11, square(exp(_1)));
  DEFINE_FUNCTOR(g22, square(exp(_1)));
  DEFINE_FUNCTOR(g33, square(exp(_1) * sin(_theta)));

  HD_INLINE void PosToCartesian(Vec3<Scalar>& pos) const {
    // detail::general_sph_pos_to_cart(pos, exp(_1), g33);
    detail::log_sph_pos_to_cart(pos);
  }

  HD_INLINE void PosFromCartesian(Vec3<Scalar>& pos) const {
    detail::log_sph_pos_from_cart(pos);
  }

  // Transform a vector to Cartesian coordinates
  HD_INLINE void VectorToCartesian(Vec3<Scalar>& v,
                                   const Vec3<Scalar>& pos) const {
    detail::general_sph_vec_to_cart(v, pos);
  }

  // Transform a vector from Cartesian coordinates
  HD_INLINE void VectorFromCartesian(Vec3<Scalar>& v,
                                     const Vec3<Scalar>& pos) const {
    detail::general_sph_vec_from_cart(v, pos);
  }

#if defined(__AVX2__) && (defined(__ICC) || defined(__INTEL_COMPILER))
  void PosToCartesian(__m256d* pos_x, __m256d* pos_y, __m256d* pos_z) const {
    detail::log_sph_pos_to_cart(pos_x, pos_y, pos_z);
  }

  void PosFromCartesian(__m256d* pos_x, __m256d* pos_y, __m256d* pos_z) const {
    detail::log_sph_pos_from_cart(pos_x, pos_y, pos_z);
  }

  void VectorToCartesian(__m256d* v1, __m256d* v2, __m256d* v3,
                         const __m256d* pos_x, const __m256d* pos_y, const __m256d* pos_z) const {
    detail::general_sph_vec_to_cart(v1, v2, v3, pos_x, pos_y, pos_z);
  }

  void VectorFromCartesian(__m256d* v1, __m256d* v2, __m256d* v3,
                           const __m256d* pos_x, const __m256d* pos_y, const __m256d* pos_z) const {
    detail::general_sph_vec_from_cart(v1, v2, v3, pos_x, pos_y, pos_z);
  }
#endif

  bool parse(const std::string& str) {
    std::stringstream ss(str);
    std::string n = "";
    ss >> n;
    std::transform(n.begin(), n.end(), n.begin(), ::tolower);
    if (n == "log_spherical") return true;
    else return false;
  }
};

extern metric_log_spherical g_metric_log_spherical;

}

}



#endif  // _METRIC_LOG_SPHERICAL_H_
