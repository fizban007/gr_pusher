#ifndef _METRIC_CARTESIAN_H_
#define _METRIC_CARTESIAN_H_

#include "metric_base.h"

namespace Aperture {

namespace metric {

struct metric_cartesian : public metric_base<metric_cartesian> {
  enum { isLinear = 1 };
  enum { swapZY = 0 };
  enum { isGR = 0 };
  enum { isDiagonal = 1 };
  enum { canMapCartesian = 1 };
  enum { type = (int)MetricType::Cartesian };

  std::string name() const { return "Cartesian"; }

  HD_INLINE void PosToCartesian(Vec3<Scalar>& pos) const {}

  HD_INLINE void PosFromCartesian(Vec3<Scalar>& pos) const {}

  HD_INLINE void VectorToCartesian(Vec3<Scalar>& v,
                                   const Vec3<Scalar>& pos) const {}

  HD_INLINE void VectorFromCartesian(Vec3<Scalar>& v,
                                     const Vec3<Scalar>& pos) const {}

#if defined(__AVX2__) && (defined(__ICC) || defined(__INTEL_COMPILER))
  void PosToCartesian(__m256d* pos_x, __m256d* pos_y, __m256d* pos_z) const {}

  void PosFromCartesian(__m256d* pos_x, __m256d* pos_y, __m256d* pos_z) const {}

  void VectorToCartesian(__m256d* v1, __m256d* v2, __m256d* v3,
                         const __m256d* pos_x, const __m256d* pos_y, const __m256d* pos_z) const {}

  void VectorFromCartesian(__m256d* v1, __m256d* v2, __m256d* v3,
                           const __m256d* pos_x, const __m256d* pos_y, const __m256d* pos_z) const {}
#endif

  bool parse(const std::string& str) {
    std::stringstream ss(str);
    std::string n = "";
    ss >> n;
    std::transform(n.begin(), n.end(), n.begin(), ::tolower);
    if (n == "cartesian") return true;
    else return false;
  }
};

extern metric_cartesian g_metric_cartesian;

}

}

#endif  // _METRIC_CARTESIAN_H_
