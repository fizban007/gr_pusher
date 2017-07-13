#ifndef _METRIC_SPHERICAL_H_
#define _METRIC_SPHERICAL_H_

#include "metric_base.h"

namespace Aperture {

namespace metric {

using CudaLE::ConstOp;
using CudaLE::placeholders::spherical::_r;
using CudaLE::placeholders::spherical::_theta;
using CudaLE::placeholders::spherical::_phi;

struct metric_spherical : public metric_base<metric_spherical> {
  enum { isLinear = 0 };
  enum { swapZY = 1 };
  enum { isGR = 0 };
  enum { isDiagonal = 1 };
  enum { canMapCartesian = 1 };
  enum { type = (int)MetricType::spherical };

  std::string name() const { return "Spherical"; }

  DEFINE_FUNCTOR(g11, ConstOp(1.0));
  DEFINE_FUNCTOR(g22, square(_r));
  DEFINE_FUNCTOR(g33, square(_r * sin(_theta)));

  HD_INLINE void PosToCartesian(Vec3<Scalar>& pos) const {
    // detail::general_sph_pos_to_cart(pos, exp(_1), g33);
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
    ss >> n;
    std::transform(n.begin(), n.end(), n.begin(), ::tolower);
    if (n == "spherical") return true;
    else return false;
  }
};

extern metric_spherical g_metric_spherical;

}

}


#endif  // _METRIC_SPHERICAL_H_
