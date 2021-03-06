#include "CudaLE.h"
#include "fadiff.h"
#include "solve.h"
#include "utils/timer.h"
#include "vec3.h"
#include "ptc.h"
#include "metrics-analytic.h"
#include "geodesic_solve.h"
// #include "CudaLE.h"
// #include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace Aperture;
using namespace CudaLE;
// using namespace fadbad;

// #define CONNECTION(METRIC, I, A, B, X)            \
//   (METRIC.inv_g##I##0(X[0], X[1], X[2]) *         \
//        (D<B>(METRIC.g##A##0)(X[0], X[1], X[2]) +  \
//         D<A>(METRIC.g##B##0)(X[0], X[1], X[2])) + \
//    METRIC.inv_g##I##1(X[0], X[1], X[2]) *         \
//        (D<B>(METRIC.g##A##1)(X[0], X[1], X[2]) +  \
//         D<A>(METRIC.g##B##1)(X[0], X[1], X[2]) -  \
//         D<1>(METRIC.g##A##B)(X[0], X[1], X[2])) + \
//    METRIC.inv_g##I##2(X[0], X[1], X[2]) *         \
//        (D<B>(METRIC.g##A##2)(X[0], X[1], X[2]) +  \
//         D<A>(METRIC.g##B##2)(X[0], X[1], X[2]) -  \
//         D<2>(METRIC.g##A##B)(X[0], X[1], X[2])) + \
//    METRIC.inv_g##I##3(X[0], X[1], X[2]) *         \
//        (D<B>(METRIC.g##A##3)(X[0], X[1], X[2]) +  \
//         D<A>(METRIC.g##B##3)(X[0], X[1], X[2])))

// #define COEF(METRIC, I, A, B, X)                     \
//   (0.5 * (D<I>(METRIC.g00)(X[0], X[1], X[2]) *       \
//               METRIC.inv_g##A##0(X[0], X[1], X[2]) * \
//               METRIC.inv_g##B##0(X[0], X[1], X[2]) + \
//           D<I>(METRIC.g01)(X[0], X[1], X[2]) *       \
//               METRIC.inv_g##A##0(X[0], X[1], X[2]) * \
//               METRIC.inv_g##B##1(X[0], X[1], X[2]) + \
//           D<I>(METRIC.g02)(X[0], X[1], X[2]) *       \
//               METRIC.inv_g##A##0(X[0], X[1], X[2]) * \
//               METRIC.inv_g##B##2(X[0], X[1], X[2]) + \
//           D<I>(METRIC.g03)(X[0], X[1], X[2]) *       \
//               METRIC.inv_g##A##0(X[0], X[1], X[2]) * \
//               METRIC.inv_g##B##3(X[0], X[1], X[2]) + \
//           D<I>(METRIC.g10)(X[0], X[1], X[2]) *       \
//               METRIC.inv_g##A##1(X[0], X[1], X[2]) * \
//               METRIC.inv_g##B##0(X[0], X[1], X[2]) + \
//           D<I>(METRIC.g11)(X[0], X[1], X[2]) *       \
//               METRIC.inv_g##A##1(X[0], X[1], X[2]) * \
//               METRIC.inv_g##B##1(X[0], X[1], X[2]) + \
//           D<I>(METRIC.g12)(X[0], X[1], X[2]) *       \
//               METRIC.inv_g##A##1(X[0], X[1], X[2]) * \
//               METRIC.inv_g##B##2(X[0], X[1], X[2]) + \
//           D<I>(METRIC.g13)(X[0], X[1], X[2]) *       \
//               METRIC.inv_g##A##1(X[0], X[1], X[2]) * \
//               METRIC.inv_g##B##3(X[0], X[1], X[2]) + \
//           D<I>(METRIC.g20)(X[0], X[1], X[2]) *       \
//               METRIC.inv_g##A##2(X[0], X[1], X[2]) * \
//               METRIC.inv_g##B##0(X[0], X[1], X[2]) + \
//           D<I>(METRIC.g21)(X[0], X[1], X[2]) *       \
//               METRIC.inv_g##A##2(X[0], X[1], X[2]) * \
//               METRIC.inv_g##B##1(X[0], X[1], X[2]) + \
//           D<I>(METRIC.g22)(X[0], X[1], X[2]) *       \
//               METRIC.inv_g##A##2(X[0], X[1], X[2]) * \
//               METRIC.inv_g##B##2(X[0], X[1], X[2]) + \
//           D<I>(METRIC.g23)(X[0], X[1], X[2]) *       \
//               METRIC.inv_g##A##2(X[0], X[1], X[2]) * \
//               METRIC.inv_g##B##3(X[0], X[1], X[2]) + \
//           D<I>(METRIC.g30)(X[0], X[1], X[2]) *       \
//               METRIC.inv_g##A##3(X[0], X[1], X[2]) * \
//               METRIC.inv_g##B##0(X[0], X[1], X[2]) + \
//           D<I>(METRIC.g31)(X[0], X[1], X[2]) *       \
//               METRIC.inv_g##A##3(X[0], X[1], X[2]) * \
//               METRIC.inv_g##B##1(X[0], X[1], X[2]) + \
//           D<I>(METRIC.g32)(X[0], X[1], X[2]) *       \
//               METRIC.inv_g##A##3(X[0], X[1], X[2]) * \
//               METRIC.inv_g##B##2(X[0], X[1], X[2]) + \
//           D<I>(METRIC.g33)(X[0], X[1], X[2]) *       \
//               METRIC.inv_g##A##3(X[0], X[1], X[2]) * \
//               METRIC.inv_g##B##3(X[0], X[1], X[2])))

// typedef fadbad::F<double, 6> var;

// template <typename Data>
// Vec3<var>
// mid_point(const Vec3<var>& x, const Vec3<var>& x0) {
//   // TODO: check the correctness of this function
//   Vec3<var> result{0.0, 0.0, 0.0};
//   for (int i = 0; i < 3; i++) result[i] = 0.5 * (x[i] + x0[i]);
//   return result;
// }

// // template <typename Data>
// var
// quadratic_solve(const var& a, const var& b, const var& c) {
//   return (-b - sqrt(b * b - 4.0 * a * c)) / (2.0 * a);
// }

// template <typename Metric>
// var
// Gamma(const Vec3<var>& x, const Vec3<var>& u,
//       const Metric& metric) {
//   var result = 1.0;
//   // TODO: Make this more general
//   result += u[0] * u[0] / metric.g11(x[0], x[1], x[2]);
//   result += u[1] * u[1] / metric.g22(x[0], x[1], x[2]);
//   result += u[2] * u[2] / metric.g33(x[0], x[1], x[2]);
//   result = sqrt(result) / metric.alpha(x[0], x[1], x[2]);
//   // Data a = metric.g00(x[0], x[1], x[2]);
//   // Data b = 2.0 * metric.g03(x[0], x[1], x[2]) * u[2];
//   // Data c = 1.0 + metric.g11(x[0], x[1], x[2]) * u[0] * u[0] +
//   //          metric.g22(x[0], x[1], x[2]) * u[1] * u[1] +
//   //          metric.g33(x[0], x[1], x[2]) * u[2] * u[2];
//   // return quadratic_solve(a, b, c);
//   return result;
// }

// // template <typename Data>
// template <typename Metric>
// var
// u0_energy(const Vec3<var>& x, const Vec3<var>& u,
//           const Metric& metric) {
//   // return -metric.alpha(x[0], x[1], x[2]) * metric.alpha(x[0], x[1], x[2]) *
//   //        Gamma(x, u, metric);
//   var gamma = Gamma(x, u, metric);
//   var result =
//       metric.b3(x[0], x[1], x[2]) * u[2] - metric.a2(x[0], x[1], x[2]) * gamma;
//   return result;
// }

// // template <typename Data>
// // Data
// // connection(int i, int a, int b, const Vec3<Data>& x,
// //            const Boyer<Data>& metric) {
// //   Data result;
// //   if (i == 1) {

// //   }
// //   return result;
// // }
// // template <typename Data>
// // Data
// // Fu(int n, const Vec3<double>& u0, const Vec3<double>& x0, const Vec3<Data>& u,
// //    const Vec3<Data>& x, const Boyer_Lindquist<Data>& metric, double dt) {
// // }

// template <typename Metric>
// var
// Fu(int n, const Vec3<var>& u0, const Vec3<var>& x0, const Vec3<var>& u,
//    const Vec3<var>& x, const Metric& metric, double dt) {
//   Vec3<var> mid_x = mid_point(x, x0);
//   Vec3<var> mid_u = mid_point(u, u0);
  // Data gamma = Gamma(mid_x, mid_u, mid_cell, grid);
  // Data u_0 = 0.5 * u0_energy(x0, u0, metric) + 0.5 * u0_energy(x, u,
  // metric);
  // std::cout << "u0 is " << u_0.x() << std::endl;
  // var result = 0.0;
  // var gamma = Gamma(mid_x, mid_u, metric);
  // var u_0 = u0_energy(mid_x, mid_u, metric);
  // var u3 = metric.inv_g30(mid_x[0], mid_x[1], mid_x[2]) * u_0 +
  //           metric.inv_g33(mid_x[0], mid_x[1], mid_x[2]) * mid_u[2];
  // std::cout << "Gamma is " << gamma.x() << std::endl;
  // std::cout << gamma<< std::endl;

//   if (n == 0) {
//     // result += COEF(metric, 1, 0, 0, mid_x) * u_0 * u_0 / gamma +
//     //           2.0 * COEF(metric, 1, 0, 1, mid_x) * u_0 * mid_u[0] / gamma +
//     //           2.0 * COEF(metric, 1, 0, 2, mid_x) * u_0 * mid_u[1] / gamma +
//     //           2.0 * COEF(metric, 1, 0, 3, mid_x) * u_0 * mid_u[2] / gamma +
//     //           COEF(metric, 1, 1, 1, mid_x) * mid_u[0] * mid_u[0] / gamma +
//     //           COEF(metric, 1, 2, 2, mid_x) * mid_u[1] * mid_u[1] / gamma +
//     //           COEF(metric, 1, 3, 3, mid_x) * mid_u[2] * mid_u[2] / gamma +
//     //           2.0 * COEF(metric, 1, 1, 2, mid_x) * mid_u[0] * mid_u[1] /
//     //           gamma +
//     //           2.0 * COEF(metric, 1, 1, 3, mid_x) * mid_u[0] * mid_u[2] /
//     //           gamma +
//     //           2.0 * COEF(metric, 1, 2, 3, mid_x) * mid_u[1] * mid_u[2] /
//     //           gamma;
//     result += 0.5 * D<1>(metric.g00)(mid_x[0], mid_x[1], mid_x[2]) * gamma;
//     result += D<1>(metric.g03)(mid_x[0], mid_x[1], mid_x[2]) * u3;
//     result += 0.5 * D<1>(metric.g11)(mid_x[0], mid_x[1], mid_x[2]) * mid_u[0] *
//               mid_u[0] / (gamma * metric.g11(mid_x[0], mid_x[1], mid_x[2]) *
//                           metric.g11(mid_x[0], mid_x[1], mid_x[2]));
//     result += 0.5 * D<1>(metric.g22)(mid_x[0], mid_x[1], mid_x[2]) * mid_u[1] *
//               mid_u[1] / (gamma * metric.g22(mid_x[0], mid_x[1], mid_x[2]) *
//                           metric.g22(mid_x[0], mid_x[1], mid_x[2]));
//     result +=
//         0.5 * D<1>(metric.g33)(mid_x[0], mid_x[1], mid_x[2]) * u3 * u3 / gamma;
//     // result -= 0.5 * CONNECTION(metric, 1, 0, 0, mid_x) * gamma +
//     //           CONNECTION(metric, 1, 0, 1, mid_x) * mid_u[0] +
//     //           CONNECTION(metric, 1, 0, 2, mid_x) * mid_u[1] +
//     //           CONNECTION(metric, 1, 0, 3, mid_x) * mid_u[2] +
//     //           0.5 * CONNECTION(metric, 1, 1, 1, mid_x) * mid_u[0] * mid_u[0]
//     //           / gamma +
//     //           0.5 * CONNECTION(metric, 1, 2, 2, mid_x) * mid_u[1] * mid_u[1]
//     //           / gamma +
//     //           0.5 * CONNECTION(metric, 1, 3, 3, mid_x) * mid_u[2] * mid_u[2]
//     //           / gamma +
//     //           CONNECTION(metric, 1, 1, 2, mid_x) * mid_u[0] * mid_u[1] /
//     //           gamma +
//     //           CONNECTION(metric, 1, 1, 3, mid_x) * mid_u[0] * mid_u[2] /
//     //           gamma +
//     //           CONNECTION(metric, 1, 2, 3, mid_x) * mid_u[1] * mid_u[2] /
//     //           gamma;
//     // std::cout << "Result is " << result.x() << std::endl;
//     // std::cout << "Gamma is " << gamma.x() << std::endl;
//     // std::cout << (CONNECTION(metric, 1, 0, 0, mid_x)).x() << std::endl;
//     // std::cout << (CONNECTION(metric, 1, 0, 1, mid_x)).x() << std::endl;
//     // std::cout << (CONNECTION(metric, 1, 0, 2, mid_x)).x() << std::endl;
//     // std::cout << (CONNECTION(metric, 1, 0, 3, mid_x)).x() << std::endl;
//     // std::cout << (CONNECTION(metric, 1, 1, 1, mid_x)).x() << std::endl;
//     // std::cout << (CONNECTION(metric, 1, 2, 2, mid_x) * mid_u[1] * mid_u[1] /
//     // gamma).x() << std::endl;
//     // std::cout << (CONNECTION(metric, 1, 3, 3, mid_x) * mid_u[2] * mid_u[2] /
//     // gamma).x() << std::endl;
//     // std::cout << (CONNECTION(metric, 1, 1, 2, mid_x)).x() << std::endl;
//     // std::cout << (CONNECTION(metric, 1, 1, 3, mid_x)).x() << std::endl;
//     // std::cout << (CONNECTION(metric, 1, 3, 2, mid_x)).x() << std::endl;
//   } else if (n == 1) {
//     result += 0.5 * D<2>(metric.g00)(mid_x[0], mid_x[1], mid_x[2]) * gamma;
//     result += D<2>(metric.g03)(mid_x[0], mid_x[1], mid_x[2]) * u3;
//     result += 0.5 * D<2>(metric.g11)(mid_x[0], mid_x[1], mid_x[2]) * mid_u[0] *
//               mid_u[0] / (gamma * metric.g11(mid_x[0], mid_x[1], mid_x[2]) *
//                           metric.g11(mid_x[0], mid_x[1], mid_x[2]));
//     result += 0.5 * D<2>(metric.g22)(mid_x[0], mid_x[1], mid_x[2]) * mid_u[1] *
//               mid_u[1] / (gamma * metric.g22(mid_x[0], mid_x[1], mid_x[2]) *
//                           metric.g22(mid_x[0], mid_x[1], mid_x[2]));
//     result +=
//         0.5 * D<2>(metric.g33)(mid_x[0], mid_x[1], mid_x[2]) * u3 * u3 / gamma;
//     // result += COEF(metric, 2, 0, 0, mid_x) * u_0 * u_0 / gamma +
//     //           2.0 * COEF(metric, 2, 0, 1, mid_x) * u_0 * mid_u[0] / gamma +
//     //           2.0 * COEF(metric, 2, 0, 2, mid_x) * u_0 * mid_u[1] / gamma +
//     //           2.0 * COEF(metric, 2, 0, 3, mid_x) * u_0 * mid_u[2] / gamma +
//     //           COEF(metric, 2, 1, 1, mid_x) * mid_u[0] * mid_u[0] / gamma +
//     //           COEF(metric, 2, 2, 2, mid_x) * mid_u[1] * mid_u[1] / gamma +
//     //           COEF(metric, 2, 3, 3, mid_x) * mid_u[2] * mid_u[2] / gamma +
//     //           2.0 * COEF(metric, 2, 1, 2, mid_x) * mid_u[0] * mid_u[1] /
//     //           gamma +
//     //           2.0 * COEF(metric, 2, 1, 3, mid_x) * mid_u[0] * mid_u[2] /
//     //           gamma +
//     //           2.0 * COEF(metric, 2, 2, 3, mid_x) * mid_u[1] * mid_u[2] /
//     //           gamma;
//     // result -= 0.5 * CONNECTION(metric, 2, 0, 0, mid_x) * gamma +
//     //           CONNECTION(metric, 2, 0, 1, mid_x) * mid_u[0] +
//     //           CONNECTION(metric, 2, 0, 2, mid_x) * mid_u[1] +
//     //           CONNECTION(metric, 2, 0, 3, mid_x) * mid_u[2] +
//     //           0.5 * CONNECTION(metric, 2, 1, 1, mid_x) * mid_u[0] * mid_u[0]
//     //           / gamma +
//     //           0.5 * CONNECTION(metric, 2, 2, 2, mid_x) * mid_u[1] * mid_u[1]
//     //           / gamma +
//     //           0.5 * CONNECTION(metric, 2, 3, 3, mid_x) * mid_u[2] * mid_u[2]
//     //           / gamma +
//     //           CONNECTION(metric, 2, 1, 2, mid_x) * mid_u[0] * mid_u[1] /
//     //           gamma +
//     //           CONNECTION(metric, 2, 1, 3, mid_x) * mid_u[0] * mid_u[2] /
//     //           gamma +
//     //           CONNECTION(metric, 2, 2, 3, mid_x) * mid_u[1] * mid_u[2] /
//     //           gamma;
//   } else {
//     // result += COEF(metric, 3, 0, 0, mid_x) * u_0 * u_0 / gamma +
//     //           2.0 * COEF(metric, 3, 0, 1, mid_x) * u_0 * mid_u[0] / gamma +
//     //           2.0 * COEF(metric, 3, 0, 2, mid_x) * u_0 * mid_u[1] / gamma +
//     //           2.0 * COEF(metric, 3, 0, 3, mid_x) * u_0 * mid_u[2] / gamma +
//     //           COEF(metric, 3, 1, 1, mid_x) * mid_u[0] * mid_u[0] / gamma +
//     //           COEF(metric, 3, 2, 2, mid_x) * mid_u[1] * mid_u[1] / gamma +
//     //           COEF(metric, 3, 3, 3, mid_x) * mid_u[2] * mid_u[2] / gamma +
//     //           2.0 * COEF(metric, 3, 1, 2, mid_x) * mid_u[0] * mid_u[1] /
//     //           gamma +
//     //           2.0 * COEF(metric, 3, 1, 3, mid_x) * mid_u[0] * mid_u[2] /
//     //           gamma +
//     //           2.0 * COEF(metric, 3, 2, 3, mid_x) * mid_u[1] * mid_u[2] /
//     //           gamma;
//     result += 0.0;
//     // result -= 0.5 * CONNECTION(metric, 3, 0, 0, mid_x) * gamma +
//     //           CONNECTION(metric, 3, 0, 1, mid_x) * mid_u[0] +
//     //           CONNECTION(metric, 3, 0, 2, mid_x) * mid_u[1] +
//     //           CONNECTION(metric, 3, 0, 3, mid_x) * mid_u[2] +
//     //           0.5 * CONNECTION(metric, 3, 1, 1, mid_x) * mid_u[0] * mid_u[0]
//     //           / gamma +
//     //           0.5 * CONNECTION(metric, 3, 2, 2, mid_x) * mid_u[1] * mid_u[1]
//     //           / gamma +
//     //           0.5 * CONNECTION(metric, 3, 3, 3, mid_x) * mid_u[2] * mid_u[2]
//     //           / gamma +
//     //           CONNECTION(metric, 3, 1, 2, mid_x) * mid_u[0] * mid_u[1] /
//     //           gamma +
//     //           CONNECTION(metric, 3, 1, 3, mid_x) * mid_u[0] * mid_u[2] /
//     //           gamma +
//     //           CONNECTION(metric, 3, 2, 3, mid_x) * mid_u[1] * mid_u[2] /
//     //           gamma;
//   }
//   return u[n] - u0[n] - result * dt;
// }

// template <typename Metric>
// var
// Fx(int n, const Vec3<var>& u0, const Vec3<var>& x0, const Vec3<var>& u,
//    const Vec3<var>& x, const Metric& metric, double dt) {
//   Vec3<var> mid_x = mid_point(x, x0);
//   Vec3<var> mid_u = mid_point(u, u0);

//   var gamma = Gamma(mid_x, mid_u, metric);
//   if (n == 0) {
//     return x[0] - x0[0] -
//            dt * mid_u[0] / (gamma * metric.g11(mid_x[0], mid_x[1], mid_x[2]));
//   } else if (n == 1) {
//     return x[1] - x0[1] -
//            dt * mid_u[1] / (gamma * metric.g22(mid_x[0], mid_x[1], mid_x[2]));
//     // dt * mid_u[1] / gamma;
//   } else if (n == 2) {
//     return x[2] - x0[2] -
//            dt * mid_u[2] / (gamma * metric.g33(mid_x[0], mid_x[1], mid_x[2])) +
//            dt * metric.b3(mid_x[0], mid_x[1], mid_x[2]);
//     // dt * mid_u[2] / gamma;
//   } else {
//     return 0.0;
//   }
// }

// template <typename Metric>
// int
// iterate_newton(Particle<var>& p, const Metric& metric, double dt) {
//   const int max_steps = 100;
//   const double tolerance = 1.0e-10;

//   Particle<var> p0 = p;

//   for (int i = 0; i < max_steps; i++) {
//     // Initialize guess solution
//     Vec3<var> x(p.x[0], p.x[1], p.x[2]);
//     Vec3<var> u(p.u[0], p.u[1], p.u[2]);
//     x[0].diff(0);
//     x[1].diff(1);
//     x[2].diff(2);
//     u[0].diff(3);
//     u[1].diff(4);
//     u[2].diff(5);
//     std::array<var, 6> f;
//     // Compute the residual function
//     f[0] = Fx(0, p0.u, p0.x, u, x, metric, dt);
//     f[1] = Fx(1, p0.u, p0.x, u, x, metric, dt);
//     f[2] = Fx(2, p0.u, p0.x, u, x, metric, dt);
//     f[3] = Fu(0, p0.u, p0.x, u, x, metric, dt);
//     f[4] = Fu(1, p0.u, p0.x, u, x, metric, dt);
//     f[5] = Fu(2, p0.u, p0.x, u, x, metric, dt);
//     // Fill the Jacobi matrix
//     std::array<std::array<double, 6>, 6> J;
//     std::array<double, 6> F, new_f;
//     for (int m = 0; m < 6; m++) {
//       for (int n = 0; n < 6; n++) {
//         J[n][m] = f[n].d(m);
//       }
//       F[m] = f[m].x();
//     }
//     // std::cout << "(" << F[0] << ", " << F[1] << ", " << F[2] << ", " << F[3]
//     //           << ", " << F[4] << ", " << F[5] << ")" << std::endl;
//     solve(J, F, new_f);
//     p.x[0] -= new_f[0];
//     p.x[1] -= new_f[1];
//     p.x[2] -= new_f[2];
//     p.u[0] -= new_f[3];
//     p.u[1] -= new_f[4];
//     p.u[2] -= new_f[5];

//     if (norm(new_f) <= tolerance) break;
//     if (i == max_steps - 1) return 1;
//   }
//   return 0;
// }

int
main(int argc, char* argv[]) {
  SolverType type = SolverType::implicit;
  Particle<var> ptc;
  double m = 1.0;
  double a = 0.8;
  Boyer_Lindquist<var> metric(m, a);
  // Boyer_Lindquist<double> metric_num(m, a);

  // ptc.x[1] = 3.1415926535898 * 0.5;
  // double L = sqrt(3.0) * 2.0 * m + 1.0e-12;
  // ptc.u[0] = 0.0;
  // double r = L * L / (2.0 * m) + sqrt(L * L / (4.0 * m * m) - 3.0) * L;
  double r = 3.0;
  double theta = 3.1415926535898 * 0.5;
  double L = sqrt(m * r) - 2.0 * a * m / r + sqrt(m / (r * r * r)) * a * a;
  L /= sqrt(1.0 - 3.0 * m / r + 2.0 * a * sqrt(m / (r * r * r)));
  // L += 0.5;
  // double u0 = -1.0 + 2.0 / r - L * L * (4.0 / (r * r) - 8.0 / (r * r * r));
  ptc.x[0] = r;
  ptc.x[1] = theta;
  // ptc.u[2] = metric_num.inv_g30(r, theta, 0.0) * u0 + metric_num.inv_g33(r,
  // theta, 0.0) * L;
  ptc.u[2] = L;
  // ptc.u[1] = 1.0e-12 * metric_num.g11(ptc.x[0], ptc.x[1], ptc.x[2]);
  ptc.u[1] = 0.3;
  auto u_0 = u0_energy(ptc.x, ptc.u, metric, ptc.is_photon);
  // auto pos = ptc.x;
  std::cout << "Initial: " << ptc << ", u_0 = " << u_0.x() << std::endl;
  // std::cout << "Initial u0 is " << u_0 << std::endl;
  // ptc.x[1] = 10.0;

  const double dt = 0.1;
  // Particle ptc_initial = ptc;
  std::ofstream out("boyer-isco-0.1.txt", std::ofstream::out);
  out << std::setprecision(std::numeric_limits<double>::digits10 + 2);
  timer::stamp();

  int N = 10000;
  for (int n = 0; n < N; n++) {
    // std::cout << "At timestep " << n << std::endl;

    // double u_0 = u0_energy(ptc.x, ptc.u, metric_num);
    // auto pos = grid.mesh().pos_particle(p.cell, p.x);
    // auto pos = ptc.x;
    // std::cout << pos << " " << ptc.u << " " << u_0 << std::endl;

    //   out << pos.x << ", " << pos.y << ", " << pos.z << ", " << u_0
    //       << std::endl;
    int iter;
    if ((iter = iterate_newton(ptc, metric, dt, type)) == -1) {
      std::cout << "Iteration reached end without converging!" << std::endl;
      break;
    }
    if (n % 10 == 0)
      std::cout << ptc << " in " << iter << " iterations" << std::endl;
  }
  auto u_1 = u0_energy(ptc.x, ptc.u, metric, ptc.is_photon);
  // auto pos = grid.mesh().pos_particle(p.cell, p.x);
  // pos = ptc.x;
  std::cout << "Final: " << ptc << ", u_0 = " << u_1.x()
            << ", \\Delta u_0/u_0 = " << ((u_1 - u_0) / u_0).x() << std::endl;
  std::cout << "Evolved for " << N << " timesteps." << std::endl;
  timer::show_duration_since_stamp("Evolution time", "ms");
  out.close();

  return 0;
}
