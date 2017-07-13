#ifndef _METRIC_SCHWARZSCHILD_H_
#define _METRIC_SCHWARZSCHILD_H_

#include "metric_base.h"

namespace Aperture {

namespace metric {

using CudaLE::placeholders::spherical::_r;
using CudaLE::placeholders::spherical::_theta;
using CudaLE::placeholders::spherical::_phi;

struct metric_schwarzschild : public metric_base<metric_schwarzschild> {
 public:
  enum { isLinear = 0 };
  enum { swapZY = 1 };
  enum { isGR = 1 };
  enum { isDiagonal = 1 };
  enum { canMapCartesian = 0 };
  enum { type = (int)MetricType::Schwarzschild };

  double m_rg;  // r_g is the Schwarzschild radius of the metric

  std::string name() const { return "Schwarzschild" + std::to_string(m_rg); }

  DEFINE_FUNCTOR(g11, 1.0 / (1.0 - m_rg / _r));
  DEFINE_FUNCTOR(g22, _r * _r);
  DEFINE_FUNCTOR(g33, _r * _r * square(sin(_theta)));

  DEFINE_FUNCTOR(a, sqrt(1.0 - m_rg / _r));
  DEFINE_FUNCTOR(g00, CudaLE::ConstOp(-1.0) * (1.0 - m_rg / _r));

  metric_schwarzschild() : m_rg(2.0) {}
  explicit metric_schwarzschild(double rg) : m_rg(rg) {}
  ~metric_schwarzschild() {}

  HD_INLINE void PosToCartesian(Vec3<Scalar>& pos) const {
    // detail::general_sph_pos_to_cart(pos, _1, g33);
    detail::general_sph_pos_to_cart(pos);
  }

  HD_INLINE void PosFromCartesian(Vec3<Scalar>& pos) const {
    detail::general_sph_pos_from_cart(pos);
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

  bool parse(const std::string& str) {
    std::stringstream ss(str);
    std::string n = "";
    double rg = 2.0;
    ss >> n >> rg;
    std::transform(n.begin(), n.end(), n.begin(), ::tolower);
    if (n == "schwarzschild") {
      if (rg > 0.0)
        m_rg = rg;
      return true;
    }
    else return false;
  }
};

extern metric_schwarzschild g_metric_schwarzschild;

}

}


#endif  // _METRIC_SCHWARZSCHILD_H_
