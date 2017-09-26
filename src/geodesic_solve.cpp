#include "geodesic_solve.h"
#include "CudaLE.h"
#include "metrics-analytic.h"
#include "solve.h"
#include <array>

using namespace CudaLE;

#define INSTANTIATE_TEMPLATES(METRIC)                                          \
  template var Gamma<METRIC>(const Vec3<var>& x, const Vec3<var>& u,           \
                             const METRIC& metric);                            \
  template var u0_energy<METRIC>(const Vec3<var>& x, const Vec3<var>& u,       \
                                 const METRIC& metric);                        \
  template var FuncX::operator()<METRIC>(                                      \
      int n, const Vec3<var>& x0, const Vec3<var>& u0, const Vec3<var>& x,     \
      const Vec3<var>& u, const METRIC& metric, double dt);                    \
  template var FuncU::operator()<METRIC>(                                      \
      int n, const Vec3<var>& x0, const Vec3<var>& u0, const Vec3<var>& x,     \
      const Vec3<var>& u, const METRIC& metric, double dt);                    \
  template int iterate_newton<METRIC>(Particle<var> & p, const METRIC& metric, \
                                      double dt)
#define CONN(METRIC, I, A, B, X)                                               \
  D<I>(METRIC.inv_g##A##B + METRIC.b##A * METRIC.b##B / METRIC.a2)(X[0], X[1], \
                                                                   X[2])

#define COEF(METRIC, I, A, B, X)                     \
  (0.5 * (D<I>(METRIC.g00)(X[0], X[1], X[2]) *       \
              METRIC.inv_g##A##0(X[0], X[1], X[2]) * \
              METRIC.inv_g##B##0(X[0], X[1], X[2]) + \
          D<I>(METRIC.g01)(X[0], X[1], X[2]) *       \
              METRIC.inv_g##A##0(X[0], X[1], X[2]) * \
              METRIC.inv_g##B##1(X[0], X[1], X[2]) + \
          D<I>(METRIC.g02)(X[0], X[1], X[2]) *       \
              METRIC.inv_g##A##0(X[0], X[1], X[2]) * \
              METRIC.inv_g##B##2(X[0], X[1], X[2]) + \
          D<I>(METRIC.g03)(X[0], X[1], X[2]) *       \
              METRIC.inv_g##A##0(X[0], X[1], X[2]) * \
              METRIC.inv_g##B##3(X[0], X[1], X[2]) + \
          D<I>(METRIC.g10)(X[0], X[1], X[2]) *       \
              METRIC.inv_g##A##1(X[0], X[1], X[2]) * \
              METRIC.inv_g##B##0(X[0], X[1], X[2]) + \
          D<I>(METRIC.g11)(X[0], X[1], X[2]) *       \
              METRIC.inv_g##A##1(X[0], X[1], X[2]) * \
              METRIC.inv_g##B##1(X[0], X[1], X[2]) + \
          D<I>(METRIC.g12)(X[0], X[1], X[2]) *       \
              METRIC.inv_g##A##1(X[0], X[1], X[2]) * \
              METRIC.inv_g##B##2(X[0], X[1], X[2]) + \
          D<I>(METRIC.g13)(X[0], X[1], X[2]) *       \
              METRIC.inv_g##A##1(X[0], X[1], X[2]) * \
              METRIC.inv_g##B##3(X[0], X[1], X[2]) + \
          D<I>(METRIC.g20)(X[0], X[1], X[2]) *       \
              METRIC.inv_g##A##2(X[0], X[1], X[2]) * \
              METRIC.inv_g##B##0(X[0], X[1], X[2]) + \
          D<I>(METRIC.g21)(X[0], X[1], X[2]) *       \
              METRIC.inv_g##A##2(X[0], X[1], X[2]) * \
              METRIC.inv_g##B##1(X[0], X[1], X[2]) + \
          D<I>(METRIC.g22)(X[0], X[1], X[2]) *       \
              METRIC.inv_g##A##2(X[0], X[1], X[2]) * \
              METRIC.inv_g##B##2(X[0], X[1], X[2]) + \
          D<I>(METRIC.g23)(X[0], X[1], X[2]) *       \
              METRIC.inv_g##A##2(X[0], X[1], X[2]) * \
              METRIC.inv_g##B##3(X[0], X[1], X[2]) + \
          D<I>(METRIC.g30)(X[0], X[1], X[2]) *       \
              METRIC.inv_g##A##3(X[0], X[1], X[2]) * \
              METRIC.inv_g##B##0(X[0], X[1], X[2]) + \
          D<I>(METRIC.g31)(X[0], X[1], X[2]) *       \
              METRIC.inv_g##A##3(X[0], X[1], X[2]) * \
              METRIC.inv_g##B##1(X[0], X[1], X[2]) + \
          D<I>(METRIC.g32)(X[0], X[1], X[2]) *       \
              METRIC.inv_g##A##3(X[0], X[1], X[2]) * \
              METRIC.inv_g##B##2(X[0], X[1], X[2]) + \
          D<I>(METRIC.g33)(X[0], X[1], X[2]) *       \
              METRIC.inv_g##A##3(X[0], X[1], X[2]) * \
              METRIC.inv_g##B##3(X[0], X[1], X[2])))

namespace Aperture {

const int max_iter = 100;
const double tolerance = 1.0e-10;

Vec3<var>
mid_point(const Vec3<var>& x, const Vec3<var>& x0) {
  Vec3<var> result;
  for (int i = 0; i < 3; i++) result[i] = 0.5 * (x[i] + x0[i]);
  return result;
}

template <typename Metric>
var
Gamma(const Vec3<var>& x, const Vec3<var>& u, const Metric& metric) {
  var result = 1.0;

  result += u[0] * u[0] * metric.inv_g11(x[0], x[1], x[2]) +
            u[1] * u[1] * metric.inv_g22(x[0], x[1], x[2]) +
            u[2] * u[2] * metric.inv_g33(x[0], x[1], x[2]) +
            2.0 * u[0] * u[1] * metric.inv_g12(x[0], x[1], x[2]) +
            2.0 * u[0] * u[2] * metric.inv_g13(x[0], x[1], x[2]) +
            2.0 * u[1] * u[2] * metric.inv_g23(x[0], x[1], x[2]);
  var ub = u[0] * metric.b1(x[0], x[1], x[2]) +
           u[1] * metric.b2(x[0], x[1], x[2]) +
           u[2] * metric.b3(x[0], x[1], x[2]);
  result += ub * ub / metric.a2(x[0], x[1], x[2]);
  return sqrt(result) / metric.alpha(x[0], x[1], x[2]);

  // return result;
}

template <typename Metric>
var
u0_energy(const Vec3<var>& x, const Vec3<var>& u, const Metric& metric) {
  var result = -Gamma(x, u, metric) * metric.a2(x[0], x[1], x[2]);
  result += metric.b1(x[0], x[1], x[2]) * u[0] +
            metric.b2(x[0], x[1], x[2]) * u[1] +
            metric.b3(x[0], x[1], x[2]) * u[2];

  return result;
}

template <typename Metric>
var
FuncX::operator()(int n, const Vec3<var>& x0, const Vec3<var>& u0,
                  const Vec3<var>& x, const Vec3<var>& u, const Metric& metric,
                  double dt) {
  // var result;
  Vec3<var> mid_x = mid_point(x, x0);
  Vec3<var> mid_u = mid_point(u, u0);

  var gamma = Gamma(mid_x, mid_u, metric);
  var ub = (mid_u[0] * metric.b1(mid_x[0], mid_x[1], mid_x[2]) +
           mid_u[1] * metric.b2(mid_x[0], mid_x[1], mid_x[2]) +
            mid_u[2] * metric.b3(mid_x[0], mid_x[1], mid_x[2])) / gamma;
  if (n == 0) {
    return x[0] - x0[0] -
           dt * ((mid_u[0] * metric.inv_g11(mid_x[0], mid_x[1], mid_x[2]) +
                  mid_u[1] * metric.inv_g12(mid_x[0], mid_x[1], mid_x[2]) +
                  mid_u[2] * metric.inv_g13(mid_x[0], mid_x[1], mid_x[2])) /
                     gamma +
                 (ub / metric.a2(mid_x[0], mid_x[1], mid_x[2]) - 1.0) *
                     metric.b1(mid_x[0], mid_x[1], mid_x[2]));
  } else if (n == 1) {
    return x[1] - x0[1] -
        dt * ((mid_u[0] * metric.inv_g21(mid_x[0], mid_x[1], mid_x[2]) +
               mid_u[1] * metric.inv_g22(mid_x[0], mid_x[1], mid_x[2]) +
               mid_u[2] * metric.inv_g23(mid_x[0], mid_x[1], mid_x[2])) /
              gamma +
              (ub / metric.a2(mid_x[0], mid_x[1], mid_x[2]) - 1.0) *
              metric.b2(mid_x[0], mid_x[1], mid_x[2]));
  } else if (n == 2) {
    return x[2] - x0[2] -
        dt * ((mid_u[0] * metric.inv_g31(mid_x[0], mid_x[1], mid_x[2]) +
               mid_u[1] * metric.inv_g32(mid_x[0], mid_x[1], mid_x[2]) +
               mid_u[2] * metric.inv_g33(mid_x[0], mid_x[1], mid_x[2])) /
              gamma +
              (ub / metric.a2(mid_x[0], mid_x[1], mid_x[2]) - 1.0) *
              metric.b3(mid_x[0], mid_x[1], mid_x[2]));
  } else {
    return 0.0;
  }

  // return result;
}

template <typename Metric>
var
FuncU::operator()(int n, const Vec3<var>& x0, const Vec3<var>& u0,
                  const Vec3<var>& x, const Vec3<var>& u, const Metric& metric,
                  double dt) {
  // var result;
  Vec3<var> mid_x = mid_point(x, x0);
  Vec3<var> mid_u = mid_point(u, u0);
  // var gamma = Gamma(mid_x, mid_u, mid_cell, grid);
  // std::cout << "Gamma is " << gamma.x() << std::endl;
  // var u_0 = 0.5 * u0_energy(x0, u0, metric) + 0.5 * u0_energy(x, u,
  // metric);
  var u_0 = u0_energy(mid_x, mid_u, metric);
  // std::cout << "u0 is " << u_0.x() << std::endl;
  var result = 0.0;
  var gamma = Gamma(mid_x, mid_u, metric);

  if (n == 0) {
    // result += COEF(metric, 1, 0, 0, mid_x) * u_0 * u_0 / gamma +
    //           2.0 * COEF(metric, 1, 0, 1, mid_x) * u_0 * mid_u[0] / gamma +
    //           2.0 * COEF(metric, 1, 0, 2, mid_x) * u_0 * mid_u[1] / gamma +
    //           2.0 * COEF(metric, 1, 0, 3, mid_x) * u_0 * mid_u[2] / gamma +
    //           COEF(metric, 1, 1, 1, mid_x) * mid_u[0] * mid_u[0] / gamma +
    //           COEF(metric, 1, 2, 2, mid_x) * mid_u[1] * mid_u[1] / gamma +
    //           COEF(metric, 1, 3, 3, mid_x) * mid_u[2] * mid_u[2] / gamma +
    //           2.0 * COEF(metric, 1, 1, 2, mid_x) * mid_u[0] * mid_u[1] /
    //           gamma + 2.0 * COEF(metric, 1, 1, 3, mid_x) * mid_u[0] *
    //           mid_u[2] / gamma + 2.0 * COEF(metric, 1, 2, 3, mid_x) *
    //           mid_u[1] * mid_u[2] / gamma;
    result += D<1>(metric.b1)(mid_x[0], mid_x[1], mid_x[2]) * mid_u[0] +
              D<1>(metric.b2)(mid_x[0], mid_x[1], mid_x[2]) * mid_u[1] +
              D<1>(metric.b3)(mid_x[0], mid_x[1], mid_x[2]) * mid_u[2];
    result -= metric.alpha(mid_x[0], mid_x[1], mid_x[2]) * gamma *
              D<1>(metric.alpha)(mid_x[0], mid_x[1], mid_x[2]);
    result -= 0.5 * CONN(metric, 1, 1, 1, mid_x) * mid_u[0] * mid_u[0] / gamma +
              0.5 * CONN(metric, 1, 2, 2, mid_x) * mid_u[1] * mid_u[1] / gamma +
              0.5 * CONN(metric, 1, 3, 3, mid_x) * mid_u[2] * mid_u[2] / gamma +
              CONN(metric, 1, 1, 2, mid_x) * mid_u[0] * mid_u[1] / gamma +
              CONN(metric, 1, 1, 3, mid_x) * mid_u[0] * mid_u[2] / gamma +
              CONN(metric, 1, 2, 3, mid_x) * mid_u[1] * mid_u[2] / gamma;
  } else if (n == 1) {
    // result += COEF(metric, 2, 0, 0, mid_x) * u_0 * u_0 / gamma +
    //           2.0 * COEF(metric, 2, 0, 1, mid_x) * u_0 * mid_u[0] / gamma +
    //           2.0 * COEF(metric, 2, 0, 2, mid_x) * u_0 * mid_u[1] / gamma +
    //           2.0 * COEF(metric, 2, 0, 3, mid_x) * u_0 * mid_u[2] / gamma +
    //           COEF(metric, 2, 1, 1, mid_x) * mid_u[0] * mid_u[0] / gamma +
    //           COEF(metric, 2, 2, 2, mid_x) * mid_u[1] * mid_u[1] / gamma +
    //           COEF(metric, 2, 3, 3, mid_x) * mid_u[2] * mid_u[2] / gamma +
    //           2.0 * COEF(metric, 2, 1, 2, mid_x) * mid_u[0] * mid_u[1] /
    //           gamma + 2.0 * COEF(metric, 2, 1, 3, mid_x) * mid_u[0] *
    //           mid_u[2] / gamma + 2.0 * COEF(metric, 2, 2, 3, mid_x) *
    //           mid_u[1] * mid_u[2] / gamma;
    result += D<2>(metric.b1)(mid_x[0], mid_x[1], mid_x[2]) * mid_u[0] +
              D<2>(metric.b2)(mid_x[0], mid_x[1], mid_x[2]) * mid_u[1] +
              D<2>(metric.b3)(mid_x[0], mid_x[1], mid_x[2]) * mid_u[2];
    result -= metric.alpha(mid_x[0], mid_x[1], mid_x[2]) * gamma *
              D<2>(metric.alpha)(mid_x[0], mid_x[1], mid_x[2]);
    result -= 0.5 * CONN(metric, 2, 1, 1, mid_x) * mid_u[0] * mid_u[0] / gamma +
              0.5 * CONN(metric, 2, 2, 2, mid_x) * mid_u[1] * mid_u[1] / gamma +
              0.5 * CONN(metric, 2, 3, 3, mid_x) * mid_u[2] * mid_u[2] / gamma +
              CONN(metric, 2, 1, 2, mid_x) * mid_u[0] * mid_u[1] / gamma +
              CONN(metric, 2, 1, 3, mid_x) * mid_u[0] * mid_u[2] / gamma +
              CONN(metric, 2, 2, 3, mid_x) * mid_u[1] * mid_u[2] / gamma;
  } else {
    // result += COEF(metric, 3, 0, 0, mid_x) * u_0 * u_0 / gamma +
    //           2.0 * COEF(metric, 3, 0, 1, mid_x) * u_0 * mid_u[0] / gamma +
    //           2.0 * COEF(metric, 3, 0, 2, mid_x) * u_0 * mid_u[1] / gamma +
    //           2.0 * COEF(metric, 3, 0, 3, mid_x) * u_0 * mid_u[2] / gamma +
    //           COEF(metric, 3, 1, 1, mid_x) * mid_u[0] * mid_u[0] / gamma +
    //           COEF(metric, 3, 2, 2, mid_x) * mid_u[1] * mid_u[1] / gamma +
    //           COEF(metric, 3, 3, 3, mid_x) * mid_u[2] * mid_u[2] / gamma +
    //           2.0 * COEF(metric, 3, 1, 2, mid_x) * mid_u[0] * mid_u[1] /
    //           gamma + 2.0 * COEF(metric, 3, 1, 3, mid_x) * mid_u[0] *
    //           mid_u[2] / gamma + 2.0 * COEF(metric, 3, 2, 3, mid_x) *
    //           mid_u[1] * mid_u[2] / gamma;
    result += D<3>(metric.b1)(mid_x[0], mid_x[1], mid_x[2]) * mid_u[0] +
              D<3>(metric.b2)(mid_x[0], mid_x[1], mid_x[2]) * mid_u[1] +
              D<3>(metric.b3)(mid_x[0], mid_x[1], mid_x[2]) * mid_u[2];
    result -= metric.alpha(mid_x[0], mid_x[1], mid_x[2]) * gamma *
              D<3>(metric.alpha)(mid_x[0], mid_x[1], mid_x[2]);
    result -= 0.5 * CONN(metric, 3, 1, 1, mid_x) * mid_u[0] * mid_u[0] / gamma +
              0.5 * CONN(metric, 3, 2, 2, mid_x) * mid_u[1] * mid_u[1] / gamma +
              0.5 * CONN(metric, 3, 3, 3, mid_x) * mid_u[2] * mid_u[2] / gamma +
              CONN(metric, 3, 1, 2, mid_x) * mid_u[0] * mid_u[1] / gamma +
              CONN(metric, 3, 1, 3, mid_x) * mid_u[0] * mid_u[2] / gamma +
              CONN(metric, 3, 2, 3, mid_x) * mid_u[1] * mid_u[2] / gamma;
  }
  return u[n] - u0[n] - result * dt;
  // return result;
}

template <typename Metric>
var
HamX::operator()(int n, const Vec3<var>& x0, const Vec3<var>& u0,
                 const Vec3<var>& x, const Vec3<var>& u, const Metric& metric,
                 double dt) {
  var result = 0.0;
  if (x == 0) {
    // result -= (dt / 6.0) * (H)
  }
  return result;
}

template <typename Metric>
var
HamU::operator()(int n, const Vec3<var>& x0, const Vec3<var>& u0,
                 const Vec3<var>& x, const Vec3<var>& u, const Metric& metric,
                 double dt) {
  var result = 0.0;
  return result;
}

template <typename Metric>
int
iterate_newton(Particle<var>& p, const Metric& metric, double dt) {
  Particle<var> p0 = p;
  FuncX Fx;
  FuncU Fu;

  for (int i = 0; i < max_iter; i++) {
    // Initialize guess solution
    Vec3<var> x(p.x[0], p.x[1], p.x[2]);
    Vec3<var> u(p.u[0], p.u[1], p.u[2]);
    x[0].diff(0);
    x[1].diff(1);
    x[2].diff(2);
    u[0].diff(3);
    u[1].diff(4);
    u[2].diff(5);
    std::array<var, 6> f;
    // Compute the residual function
    f[0] = Fx(0, p0.x, p0.u, x, u, metric, dt);
    f[1] = Fx(1, p0.x, p0.u, x, u, metric, dt);
    f[2] = Fx(2, p0.x, p0.u, x, u, metric, dt);
    f[3] = Fu(0, p0.x, p0.u, x, u, metric, dt);
    f[4] = Fu(1, p0.x, p0.u, x, u, metric, dt);
    f[5] = Fu(2, p0.x, p0.u, x, u, metric, dt);
    // Fill the Jacobi matrix
    std::array<std::array<double, 6>, 6> J;
    std::array<double, 6> F, new_f;
    for (int m = 0; m < 6; m++) {
      for (int n = 0; n < 6; n++) {
        J[n][m] = f[n].d(m);
      }
      F[m] = f[m].x();
    }
    // std::cout << "(" << F[0] << ", " << F[1] << ", " << F[2] << ", " << F[3]
    //           << ", " << F[4] << ", " << F[5] << ")" << std::endl;
    solve(J, F, new_f);
    p.x[0] -= new_f[0];
    p.x[1] -= new_f[1];
    p.x[2] -= new_f[2];
    p.u[0] -= new_f[3];
    p.u[1] -= new_f[4];
    p.u[2] -= new_f[5];

    if (norm(new_f) <= tolerance) break;
    if (i == max_iter - 1) return 1;
  }
  return 0;
}

INSTANTIATE_TEMPLATES(Schwarzschild<var>);
INSTANTIATE_TEMPLATES(Boyer_Lindquist<var>);
INSTANTIATE_TEMPLATES(Cartesian_KS<var>);
INSTANTIATE_TEMPLATES(Kerr_Schild<var>);

}  // namespace Aperture
