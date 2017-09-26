#ifndef _GEODESIC_SOLVE_H_
#define _GEODESIC_SOLVE_H_

#include "fadiff.h"
#include "vec3.h"
#include "ptc.h"

namespace Aperture {

typedef fadbad::F<double, 6> var;

Vec3<var> mid_point(const Vec3<var>& x, const Vec3<var>& x0);

template <typename Metric>
var Gamma(const Vec3<var>& x, const Vec3<var>& u, const Metric& metric, bool is_photon);

template <typename Metric>
var u0_energy(const Vec3<var>& x, const Vec3<var>& u, const Metric& metric, bool is_photon);

struct FuncX {
  template <typename Metric>
  var operator()(int n, const Vec3<var>& x0, const Vec3<var>& u0, const Vec3<var>& x,
                 const Vec3<var>& u, const Metric& metric, double dt, bool is_photon);
};

struct FuncU {
  template <typename Metric>
  var operator()(int n, const Vec3<var>& x0, const Vec3<var>& u0, const Vec3<var>& x,
                 const Vec3<var>& u, const Metric& metric, double dt, bool is_photon);
};

struct HamX {
  template <typename Metric>
  var operator()(int n, const Vec3<var>& x0, const Vec3<var>& u0, const Vec3<var>& x,
                 const Vec3<var>& u, const Metric& metric, double dt);
};

struct HamU {
  template <typename Metric>
  var operator()(int n, const Vec3<var>& x0, const Vec3<var>& u0, const Vec3<var>& x,
                 const Vec3<var>& u, const Metric& metric, double dt);
};

template <typename Metric>
int iterate_newton(Particle<var>& p, const Metric& metric, double dt);


}

#endif  // _GEODESIC_SOLVE_H_
