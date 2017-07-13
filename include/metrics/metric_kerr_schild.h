#ifndef _METRIC_KERR_SCHILD_H_
#define _METRIC_KERR_SCHILD_H_

#include "metric_base.h"

namespace Aperture {

namespace metric {

struct metric_kerr_schild : public metric_base<metric_kerr_schild> {
 public:
  double m_rg, m_a;  // r_g is 2M and a is the angular momentum of the BH
  enum { isLinear = 0 };
  enum { swapZY = 1 };
  enum { isGR = 1 };
  enum { isDiagonal = 0 };
  enum { canMapCartesian = 0 };
  enum { type = (int)MetricType::Kerr_Schild };

  std::string name() const {
    return "Kerr-Schild" + std::to_string(m_rg) + std::to_string(m_a);
  }

  DEFINE_FUNCTOR(rho2, _1 * _1 + m_a * m_a * square(cos(_2)));
  DEFINE_FUNCTOR(z, _1 * m_rg / rho2);
  DEFINE_FUNCTOR(delta, _1 * _1 + m_a * m_a - m_rg * _1);
  DEFINE_FUNCTOR(g11, 1.0 + z);
  DEFINE_FUNCTOR(g22, rho2);
  DEFINE_FUNCTOR(g33, (square(_1 * _1 + m_a * m_a) -
                       square(m_a * sin(_2)) * delta) *
                 square(sin(_2)) / g22);
  DEFINE_FUNCTOR(g13, -m_a * square(sin(_2)) * g11);
  DEFINE_FUNCTOR(g31, g13);

  DEFINE_FUNCTOR(a, 1.0 / sqrt(g11));
  DEFINE_FUNCTOR(b1, z / g11);

  DEFINE_FUNCTOR(g00, z - 1.0);
  DEFINE_FUNCTOR(g01, z);
  DEFINE_FUNCTOR(g03, -1.0 * z * m_a * square(sin(_2)));
  DEFINE_FUNCTOR(g10, z);
  DEFINE_FUNCTOR(g30, -1.0 * z * m_a * square(sin(_2)));

  metric_kerr_schild() : m_rg(2.0), m_a(0.0) {}
  metric_kerr_schild(double rg, double a) : m_rg(rg), m_a(a) {}
  ~metric_kerr_schild() {}

  HD_INLINE void PosToCartesian(Vec3<Scalar>& pos) const {
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
    double rg = 2.0, a = 0.0;
    ss >> n >> rg >> a;
    std::transform(n.begin(), n.end(), n.begin(), ::tolower);
    if (n == "kerr_schild") {
      if (rg > 0.0) m_rg = rg;
      if (a > 0.0) m_a = a;
      return true;
    }
    else return false;
  }
};

extern metric_kerr_schild g_metric_kerr_schild;

}

}

#endif  // _METRIC_KERR_SCHILD_H_
