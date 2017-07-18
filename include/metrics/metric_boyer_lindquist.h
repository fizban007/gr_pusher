#ifndef _METRIC_BOYER_LINDQUIST_H_
#define _METRIC_BOYER_LINDQUIST_H_

#include "metric_base.h"

namespace Aperture {

namespace metric {

using CudaLE::placeholders::spherical::_r;
using CudaLE::placeholders::spherical::_theta;
using CudaLE::placeholders::spherical::_phi;

struct metric_boyer_lindquist : public metric_base<metric_boyer_lindquist> {
 public:
  double m_rg, m_a;  // r_g is 2M and a is the angular momentum of the BH
  enum { isLinear = 0 };
  enum { swapZY = 1 };
  enum { isGR = 1 };
  enum { isDiagonal = 1 };
  enum { canMapCartesian = 0 };
  enum { type = (int)MetricType::Boyer_Lindquist };

  std::string name() const {
    return "Boyer-Lindquist" + std::to_string(m_rg) + std::to_string(m_a);
  }

  DEFINE_FUNCTOR(rho2, _r * _r + square(m_a * cos(_theta)));
  DEFINE_FUNCTOR(z, _r * m_rg / rho2);
  DEFINE_FUNCTOR(delta, _r* _r + m_a * m_a - m_rg * _r);
  DEFINE_FUNCTOR(sigma, square(_r * _r + m_a * m_a) -
                 delta * square(m_a * sin(_theta)));
  DEFINE_FUNCTOR(g11, rho2 / delta);
  DEFINE_FUNCTOR(g22, rho2);
  DEFINE_FUNCTOR(g33, sigma * square(sin(_theta)) / rho2);

  DEFINE_FUNCTOR(a, sqrt(rho2 * delta / sigma));
  DEFINE_FUNCTOR(b3, -m_a * m_rg * _r / sigma);

  DEFINE_FUNCTOR(g00, z - 1.0);
  DEFINE_FUNCTOR(g03, -1.0 * z * m_a * square(sin(_theta)));
  DEFINE_FUNCTOR(g30, -1.0 * z * m_a * square(sin(_theta)));

  metric_boyer_lindquist() : m_rg(2.0), m_a(0.0) {}
  metric_boyer_lindquist(double rg, double a) : m_rg(rg), m_a(a) {}
  ~metric_boyer_lindquist() {}

  HD_INLINE void PosToCartesian(Vec3<Scalar>& pos) const {
    // detail::general_sph_pos_to_cart(pos, _r, g33);
    detail::general_sph_pos_to_cart(pos, sqrt(_r * _r + m_a * m_a), (_r * _r + m_a * m_a) * square(sin(_theta)));
    // detail::general_sph_pos_to_cart(pos);
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
    double rg = 2.0, a = 0.0;
    ss >> n >> rg >> a;
    std::transform(n.begin(), n.end(), n.begin(), ::tolower);
    if (n == "boyer_lindquist") {
      if (rg > 0.0) m_rg = rg;
      if (a > 0.0) m_a = a;
      return true;
    }
    else return false;
  }
};

extern metric_boyer_lindquist g_metric_boyer_lindquist;

}


}

#endif  // _METRIC_BOYER_LINDQUIST_H_
