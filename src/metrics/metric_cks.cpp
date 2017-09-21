#include "metrics/metric_cks.h"
#include "CudaLE.h"
#include "fadiff.h"
#include "solve.h"
#include <array>

using namespace CudaLE;

#define CONNECTION(METRIC, I, A, B, X)                     \
  METRIC.inv_g##I##0(X[0], X[1], X[2]) *                   \
          (2.0 * D<B>(METRIC.g##A##0)(X[0], X[1], X[2])) + \
      METRIC.inv_g##I##1(X[0], X[1], X[2]) *               \
          (2.0 * D<B>(METRIC.g##A##1)(X[0], X[1], X[2]) -  \
           D<1>(METRIC.g##A##B)(X[0], X[1], X[2])) +       \
      METRIC.inv_g##I##2(X[0], X[1], X[2]) *               \
          (2.0 * D<B>(METRIC.g##A##2)(X[0], X[1], X[2]) -  \
           D<2>(METRIC.g##A##B)(X[0], X[1], X[2])) +       \
      METRIC.inv_g##I##3(X[0], X[1], X[2]) *               \
          (2.0 * D<B>(METRIC.g##A##3)(X[0], X[1], X[2]) -  \
           D<3>(METRIC.g##A##B)(X[0], X[1], X[2]))

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

typedef fadbad::F<double, 6> var;

template struct Cartesian_KS<double>;
template struct Cartesian_KS<var>;

template <typename Data>
Vec3<Data>
mid_point(const Vec3<Data>& x, const Vec3<double>& x0) {
  // TODO: check the correctness of this function
  Vec3<Data> result{0.0, 0.0, 0.0};
  for (int i = 0; i < 3; i++) result[i] = 0.5 * (x[i] + x0[i]);
  return result;
}

template <typename Data>
Data
quadratic_solve(const Data& a, const Data& b, const Data& c, float sign) {
  return (-b + sign * sqrt(b * b - 4.0 * a * c)) / (2.0 * a);
}

template <typename Data>
Data
Gamma(const Vec3<Data>& x, const Vec3<Data>& u,
      const Cartesian_KS<Data>& metric) {
  Data result = 1.0;
  // TODO: Make this more general
  result += u[0] * u[0] * metric.inv_g11(x[0], x[1], x[2]) +
            u[0] * u[1] * metric.inv_g12(x[0], x[1], x[2]) +
            u[0] * u[2] * metric.inv_g13(x[0], x[1], x[2]) +
            u[1] * u[0] * metric.inv_g21(x[0], x[1], x[2]) +
            u[1] * u[1] * metric.inv_g22(x[0], x[1], x[2]) +
            u[1] * u[2] * metric.inv_g23(x[0], x[1], x[2]) +
            u[2] * u[0] * metric.inv_g31(x[0], x[1], x[2]) +
            u[2] * u[1] * metric.inv_g32(x[0], x[1], x[2]) +
            u[2] * u[2] * metric.inv_g33(x[0], x[1], x[2]);
  Data ub = u[0] * metric.b1(x[0], x[1], x[2]) +
            u[1] * metric.b2(x[0], x[1], x[2]) +
            u[2] * metric.b3(x[0], x[1], x[2]);
  result += ub * ub / metric.a2(x[0], x[1], x[2]);
  return sqrt(result) / metric.alpha(x[0], x[1], x[2]);
  // result += u[0] * u[0] / metric.g11(x[0], x[1], x[2]);
  // result += u[1] * u[1] / metric.g22(x[0], x[1], x[2]);
  // result += u[2] * u[2] / metric.g33(x[0], x[1], x[2]);
  // result = sqrt(result) / metric.alpha(x[0], x[1], x[2]);
  // Data a = metric.g00(x[0], x[1], x[2]);
  // Data b = 2.0 * (metric.g01(x[0], x[1], x[2]) * u[0] + metric.g02(x[0],
  // x[1], x[2]) * u[1] +
  //                 metric.g03(x[0], x[1], x[2]) * u[2]);
  // // Data c = 1.0 + metric.g11(x[0], x[1], x[2]) * u[0] * u[0] +
  // Data c = metric.g11(x[0], x[1], x[2]) * u[0] * u[0] +
  //          metric.g22(x[0], x[1], x[2]) * u[1] * u[1] +
  //          metric.g33(x[0], x[1], x[2]) * u[2] * u[2] +
  //          2.0 * metric.g12(x[0], x[1], x[2]) * u[0] * u[1] +
  //          2.0 * metric.g13(x[0], x[1], x[2]) * u[0] * u[2] +
  //          2.0 * metric.g23(x[0], x[1], x[2]) * u[1] * u[2];
  // return quadratic_solve(a, b, c);
}

// template <typename Data1, typename Data2, typename Data3, typename Data4,
          // typename Data5, typename Data6>
// var
// Gamma(const Data1& x1, const Data2& u1, const Data3& x2, const Data4& u2,
//       const Data5& x3, const Data6& u3, const Cartesian_KS<)

template <typename Data>
Data
u0_energy(const Vec3<Data>& x, const Vec3<Data>& u,
          const Cartesian_KS<Data>& metric) {
  Data result = -Gamma(x, u, metric) * metric.a2(x[0], x[1], x[2]);
  result += metric.b1(x[0], x[1], x[2]) * u[0] +
            metric.b2(x[0], x[1], x[2]) * u[1] +
            metric.b3(x[0], x[1], x[2]) * u[2];
  return result;
  // return -metric.alpha(x[0], x[1], x[2]) * metric.alpha(x[0], x[1], x[2]) *
  //        Gamma(x, u, metric);
  // Data result = Gamma(x, u, metric);
  // Data result = metric.g00(x[0], x[1], x[2]) * Gamma(x, u, metric);
  // result += metric.g01(x[0], x[1], x[2]) * u[0] + metric.g02(x[0], x[1],
  // x[2]) * u[1] +
  //           metric.g03(x[0], x[1], x[2]) * u[2];
  // Data a = metric.inv_g00(x[0], x[1], x[2]);
  // Data b = 2.0 * (metric.inv_g01(x[0], x[1], x[2]) * u[0] +
  // metric.inv_g02(x[0], x[1], x[2]) * u[1] +
  //                 metric.inv_g03(x[0], x[1], x[2]) * u[2]);
  // // Data c = 1.0 + metric.g11(x[0], x[1], x[2]) * u[0] * u[0] +
  // Data c = metric.inv_g11(x[0], x[1], x[2]) * u[0] * u[0] +
  //          metric.inv_g22(x[0], x[1], x[2]) * u[1] * u[1] +
  //          metric.inv_g33(x[0], x[1], x[2]) * u[2] * u[2] +
  //          2.0 * metric.inv_g12(x[0], x[1], x[2]) * u[0] * u[1] +
  //          2.0 * metric.inv_g13(x[0], x[1], x[2]) * u[0] * u[2] +
  //          2.0 * metric.inv_g23(x[0], x[1], x[2]) * u[1] * u[2];
  // return quadratic_solve(a, b, c, -1.0);
  // return result;
}

template <typename Data>
Data
Fu(int n, const Vec3<double>& u0, const Vec3<double>& x0, const Vec3<Data>& u,
   const Vec3<Data>& x, const Cartesian_KS<Data>& metric, double dt) {
  Vec3<Data> mid_x = mid_point(x, x0);
  Vec3<Data> mid_u = mid_point(u, u0);
  // Data gamma = Gamma(mid_x, mid_u, mid_cell, grid);
  // std::cout << "Gamma is " << gamma.x() << std::endl;
  // Data u_0 = 0.5 * u0_energy(x0, u0, metric) + 0.5 * u0_energy(x, u,
  // metric);
  Data u_0 = u0_energy(mid_x, mid_u, metric);
  // std::cout << "u0 is " << u_0.x() << std::endl;
  Data result = 0.0;
  Data gamma = Gamma(mid_x, mid_u, metric);

  if (n == 0) {
    // result -= 0.5 * CONNECTION(metric, 1, 0, 0, mid_x) * gamma +
    //           CONNECTION(metric, 1, 0, 1, mid_x) * mid_u[0] +
    //           CONNECTION(metric, 1, 0, 2, mid_x) * mid_u[1] +
    //           CONNECTION(metric, 1, 0, 3, mid_x) * mid_u[2] +
    //           0.5 * CONNECTION(metric, 1, 1, 1, mid_x) * mid_u[0] * mid_u[0]
    //           / gamma + 0.5 * CONNECTION(metric, 1, 2, 2, mid_x) * mid_u[1] *
    //           mid_u[1] / gamma + 0.5 * CONNECTION(metric, 1, 3, 3, mid_x) *
    //           mid_u[2] * mid_u[2] / gamma + CONNECTION(metric, 1, 1, 2,
    //           mid_x) * mid_u[0] * mid_u[1] / gamma + CONNECTION(metric, 1, 1,
    //           3, mid_x) * mid_u[0] * mid_u[2] / gamma + CONNECTION(metric, 1,
    //           2, 3, mid_x) * mid_u[1] * mid_u[2] / gamma;
    result += COEF(metric, 1, 0, 0, mid_x) * u_0 * u_0 / gamma +
              2.0 * COEF(metric, 1, 0, 1, mid_x) * u_0 * mid_u[0] / gamma +
              2.0 * COEF(metric, 1, 0, 2, mid_x) * u_0 * mid_u[1] / gamma +
              2.0 * COEF(metric, 1, 0, 3, mid_x) * u_0 * mid_u[2] / gamma +
              COEF(metric, 1, 1, 1, mid_x) * mid_u[0] * mid_u[0] / gamma +
              COEF(metric, 1, 2, 2, mid_x) * mid_u[1] * mid_u[1] / gamma +
              COEF(metric, 1, 3, 3, mid_x) * mid_u[2] * mid_u[2] / gamma +
              2.0 * COEF(metric, 1, 1, 2, mid_x) * mid_u[0] * mid_u[1] / gamma +
              2.0 * COEF(metric, 1, 1, 3, mid_x) * mid_u[0] * mid_u[2] / gamma +
              2.0 * COEF(metric, 1, 2, 3, mid_x) * mid_u[1] * mid_u[2] / gamma;
  } else if (n == 1) {
    // result -= 0.5 * CONNECTION(metric, 2, 0, 0, mid_x) * gamma +
    //           CONNECTION(metric, 2, 0, 1, mid_x) * mid_u[0] +
    //           CONNECTION(metric, 2, 0, 2, mid_x) * mid_u[1] +
    //           CONNECTION(metric, 2, 0, 3, mid_x) * mid_u[2] +
    //           0.5 * CONNECTION(metric, 2, 1, 1, mid_x) * mid_u[0] * mid_u[0]
    //           / gamma + 0.5 * CONNECTION(metric, 2, 2, 2, mid_x) * mid_u[1] *
    //           mid_u[1] / gamma + 0.5 * CONNECTION(metric, 2, 3, 3, mid_x) *
    //           mid_u[2] * mid_u[2] / gamma + CONNECTION(metric, 2, 1, 2,
    //           mid_x) * mid_u[0] * mid_u[1] / gamma + CONNECTION(metric, 2, 1,
    //           3, mid_x) * mid_u[0] * mid_u[2] / gamma + CONNECTION(metric, 2,
    //           2, 3, mid_x) * mid_u[1] * mid_u[2] / gamma;
    result += COEF(metric, 2, 0, 0, mid_x) * u_0 * u_0 / gamma +
              2.0 * COEF(metric, 2, 0, 1, mid_x) * u_0 * mid_u[0] / gamma +
              2.0 * COEF(metric, 2, 0, 2, mid_x) * u_0 * mid_u[1] / gamma +
              2.0 * COEF(metric, 2, 0, 3, mid_x) * u_0 * mid_u[2] / gamma +
              COEF(metric, 2, 1, 1, mid_x) * mid_u[0] * mid_u[0] / gamma +
              COEF(metric, 2, 2, 2, mid_x) * mid_u[1] * mid_u[1] / gamma +
              COEF(metric, 2, 3, 3, mid_x) * mid_u[2] * mid_u[2] / gamma +
              2.0 * COEF(metric, 2, 1, 2, mid_x) * mid_u[0] * mid_u[1] / gamma +
              2.0 * COEF(metric, 2, 1, 3, mid_x) * mid_u[0] * mid_u[2] / gamma +
              2.0 * COEF(metric, 2, 2, 3, mid_x) * mid_u[1] * mid_u[2] / gamma;
  } else {
    // result -= 0.5 * CONNECTION(metric, 3, 0, 0, mid_x) * gamma +
    //           CONNECTION(metric, 3, 0, 1, mid_x) * mid_u[0] +
    //           CONNECTION(metric, 3, 0, 2, mid_x) * mid_u[1] +
    //           CONNECTION(metric, 3, 0, 3, mid_x) * mid_u[2] +
    //           0.5 * CONNECTION(metric, 3, 1, 1, mid_x) * mid_u[0] * mid_u[0]
    //           / gamma + 0.5 * CONNECTION(metric, 3, 2, 2, mid_x) * mid_u[1] *
    //           mid_u[1] / gamma + 0.5 * CONNECTION(metric, 3, 3, 3, mid_x) *
    //           mid_u[2] * mid_u[2] / gamma + CONNECTION(metric, 3, 1, 2,
    //           mid_x) * mid_u[0] * mid_u[1] / gamma + CONNECTION(metric, 3, 1,
    //           3, mid_x) * mid_u[0] * mid_u[2] / gamma + CONNECTION(metric, 3,
    //           2, 3, mid_x) * mid_u[1] * mid_u[2] / gamma;
    result += COEF(metric, 3, 0, 0, mid_x) * u_0 * u_0 / gamma +
              2.0 * COEF(metric, 3, 0, 1, mid_x) * u_0 * mid_u[0] / gamma +
              2.0 * COEF(metric, 3, 0, 2, mid_x) * u_0 * mid_u[1] / gamma +
              2.0 * COEF(metric, 3, 0, 3, mid_x) * u_0 * mid_u[2] / gamma +
              COEF(metric, 3, 1, 1, mid_x) * mid_u[0] * mid_u[0] / gamma +
              COEF(metric, 3, 2, 2, mid_x) * mid_u[1] * mid_u[1] / gamma +
              COEF(metric, 3, 3, 3, mid_x) * mid_u[2] * mid_u[2] / gamma +
              2.0 * COEF(metric, 3, 1, 2, mid_x) * mid_u[0] * mid_u[1] / gamma +
              2.0 * COEF(metric, 3, 1, 3, mid_x) * mid_u[0] * mid_u[2] / gamma +
              2.0 * COEF(metric, 3, 2, 3, mid_x) * mid_u[1] * mid_u[2] / gamma;
  }
  return u[n] - u0[n] - result * dt;
}

template <typename Data>
Data
Fx(int n, const Vec3<double>& u0, const Vec3<double>& x0, const Vec3<Data>& u,
   const Vec3<Data>& x, const Cartesian_KS<Data>& metric, double dt) {
  Vec3<Data> mid_x = mid_point(x, x0);
  Vec3<Data> mid_u = mid_point(u, u0);

  Data gamma = Gamma(mid_x, mid_u, metric);
  if (n == 0) {
    return x[0] - x0[0] -
           dt *
               (mid_u[0] * metric.inv_g11(mid_x[0], mid_x[1], mid_x[2]) +
                mid_u[1] * metric.inv_g12(mid_x[0], mid_x[1], mid_x[2]) +
                mid_u[2] * metric.inv_g13(mid_x[0], mid_x[1], mid_x[2])) /
               gamma;
  } else if (n == 1) {
    return x[1] - x0[1] -
           // dt * mid_u[1] / gamma;
           dt *
               (mid_u[0] * metric.inv_g21(mid_x[0], mid_x[1], mid_x[2]) +
                mid_u[1] * metric.inv_g22(mid_x[0], mid_x[1], mid_x[2]) +
                mid_u[2] * metric.inv_g23(mid_x[0], mid_x[1], mid_x[2])) /
               gamma;
  } else if (n == 2) {
    return x[2] - x0[2] -
           // dt * mid_u[2] / gamma;
           dt *
               (mid_u[0] * metric.inv_g31(mid_x[0], mid_x[1], mid_x[2]) +
                mid_u[1] * metric.inv_g32(mid_x[0], mid_x[1], mid_x[2]) +
                mid_u[2] * metric.inv_g33(mid_x[0], mid_x[1], mid_x[2])) /
               gamma;
  } else {
    return 0.0;
  }
}

// template <typename Data1, typename Data2, typename Data3, typename Data4,
//           typename Data5, typename Data6>
// Data
// Fu_H(int n, const Data1& u_1, const Data2& x_1, const Data3& u_2,
//      const Data4& x_2, const Data5& u_3, const Data6& x_3,
//      const Cartesian_KS<Data>& metric, double dt) {}

template <typename Metric>
int
iterate_newton(Particle& p, const Metric& metric, double dt) {
  const int max_steps = 100;
  const double tolerance = 1.0e-10;

  Particle p0 = p;

  for (int i = 0; i < max_steps; i++) {
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
    f[0] = Fx(0, p0.u, p0.x, u, x, metric, dt);
    f[1] = Fx(1, p0.u, p0.x, u, x, metric, dt);
    f[2] = Fx(2, p0.u, p0.x, u, x, metric, dt);
    f[3] = Fu(0, p0.u, p0.x, u, x, metric, dt);
    f[4] = Fu(1, p0.u, p0.x, u, x, metric, dt);
    f[5] = Fu(2, p0.u, p0.x, u, x, metric, dt);
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
    if (i == max_steps - 1) return 1;
  }
  return 0;
}

template var quadratic_solve<var>(const var& a, const var& b, const var& c,
                                  float sign);
template var Gamma<var>(const Vec3<var>& x, const Vec3<var>& u,
                        const Cartesian_KS<var>& metric);
template var u0_energy<var>(const Vec3<var>& x, const Vec3<var>& u,
                            const Cartesian_KS<var>& metric);
template var Fu<var>(int n, const Vec3<double>& u0, const Vec3<double>& x0,
                     const Vec3<var>& u, const Vec3<var>& x,
                     const Cartesian_KS<var>& metric, double dt);
template var Fx<var>(int n, const Vec3<double>& u0, const Vec3<double>& x0,
                     const Vec3<var>& u, const Vec3<var>& x,
                     const Cartesian_KS<var>& metric, double dt);

template double quadratic_solve<double>(const double& a, const double& b,
                                        const double& c, float sign);
template double Gamma<double>(const Vec3<double>& x, const Vec3<double>& u,
                              const Cartesian_KS<double>& metric);
template double u0_energy<double>(const Vec3<double>& x, const Vec3<double>& u,
                                  const Cartesian_KS<double>& metric);
template int iterate_newton<Cartesian_KS<var>>(Particle& p,
                                               const Cartesian_KS<var>& metric,
                                               double dt);

}  // namespace Aperture
