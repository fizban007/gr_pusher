#ifndef _GEODESIC_SOLVE_H_
#define _GEODESIC_SOLVE_H_

#include "fadiff.h"
#include "ptc.h"
#include "vec3.h"

namespace Aperture {

enum class SolverType {
  implicit, hamiltonian
};

typedef fadbad::F<double, 6> var;

Vec3<var> mid_point(const Vec3<var>& x, const Vec3<var>& x0);

template <typename Metric>
var Gamma(const Vec3<var>& x, const Vec3<var>& u, const Metric& metric,
          bool is_photon);

template <typename Metric>
var u0_energy(const Vec3<var>& x, const Vec3<var>& u, const Metric& metric,
              bool is_photon);

struct FuncX {
  template <typename Metric>
  var operator()(int n, const Vec3<var>& x0, const Vec3<var>& u0,
                 const Vec3<var>& x, const Vec3<var>& u, const Metric& metric,
                 double dt, bool is_photon);
};

struct FuncU {
  template <typename Metric>
  var operator()(int n, const Vec3<var>& x0, const Vec3<var>& u0,
                 const Vec3<var>& x, const Vec3<var>& u, const Metric& metric,
                 double dt, bool is_photon);
};

struct HamX {
  template <typename Metric>
  var operator()(int n, Vec3<var>& x0, Vec3<var>& u0, Vec3<var>& x,
                 Vec3<var>& u, const Metric& metric, double dt, bool is_photon);
};

struct HamU {
  template <typename Metric>
  var operator()(int n, Vec3<var>& x0, Vec3<var>& u0, Vec3<var>& x,
                 Vec3<var>& u, const Metric& metric, double dt, bool is_photon);
};

template <typename Metric>
int iterate_newton(Particle<var>& p, const Metric& metric, double dt, SolverType type);

template <typename Metric>
int iterate_rk4(Particle<var>& p, const Metric& metric, double dt);

}  // namespace Aperture

#endif  // _GEODESIC_SOLVE_H_
