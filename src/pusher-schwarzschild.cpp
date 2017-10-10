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

// template <typename Data>
// struct Schwarzschild {
//   CudaLE::Var<1, Data> _r;
//   CudaLE::Var<2, Data> _theta;
//   CudaLE::Var<3, Data> _phi;

//   double rg = 2.0;

//   Schwarzschild() {}
//   Schwarzschild(double r_g) : rg(r_g) {}

//   DEFINE_FUNCTOR(g00, (rg / _r - 1.0));
//   DEFINE_FUNCTOR(g11, 1.0 / (1.0 - rg / _r));
//   DEFINE_FUNCTOR(g22, _r* _r);
//   DEFINE_FUNCTOR(g33, _r* _r* square(sin(_theta)));

//   DEFINE_FUNCTOR(alpha, sqrt(1.0 - rg / _r));
// };

// template <typename Data>
// struct MetricCache {
//   Data g11, g22, g33;
//   Data gamma;
// };

// template <typename Double>
// Vec3<Double>
// mid_point(const Vec3<Double>& x, const Vec3<double>& x0) {
//   // TODO: check the correctness of this function
//   Vec3<Double> result{0.0, 0.0, 0.0};
//   for (int i = 0; i < 3; i++) result[i] = 0.5 * (x[i] + x0[i]);
//   return result;
// }

// template <typename Double>
// Double
// Gamma(const Vec3<Double>& x, const Vec3<Double>& u,
//       const Schwarzschild<Double>& metric, const MetricCache<Double>& cache) {
//   Double result = 1.0;
//   // TODO: Make this more general
//   result += u[0] * u[0] / cache.g11;
//   result += u[1] * u[1] / cache.g22;
//   result += u[2] * u[2] / cache.g33;
//   result = sqrt(result) / metric.alpha(x[0], x[1], x[2]);
//   return result;
// }

// template <typename Double>
// Double
// Gamma(const Vec3<Double>& x, const Vec3<Double>& u,
//       const Schwarzschild<Double>& metric) {
//   Double result = 1.0;
//   // TODO: Make this more general
//   result += u[0] * u[0] / metric.g11(x[0], x[1], x[2]);
//   result += u[1] * u[1] / metric.g22(x[0], x[1], x[2]);
//   result += u[2] * u[2] / metric.g33(x[0], x[1], x[2]);
//   result = sqrt(result) / metric.alpha(x[0], x[1], x[2]);
//   return result;
// }

// template <typename Double>
// Double
// u0_energy(const Vec3<Double>& x, const Vec3<Double>& u,
//           const Schwarzschild<Double>& metric,
//           const MetricCache<Double>& cache) {
//   return -metric.alpha(x[0], x[1], x[2]) * metric.alpha(x[0], x[1], x[2]) *
//          cache.gamma;
// }

// template <typename Double>
// Double
// u0_energy(const Vec3<Double>& x, const Vec3<Double>& u,
//           const Schwarzschild<Double>& metric) {
//   return -metric.alpha(x[0], x[1], x[2]) * metric.alpha(x[0], x[1], x[2]) *
//          Gamma(x, u, metric);
// }

// template <typename Double>
// Double
// Fu(int n, const Vec3<double>& u0, const Vec3<double>& x0, const Vec3<Double>& u,
//    const Vec3<Double>& x, const Schwarzschild<Double>& metric,
//    const MetricCache<Double>& cache, double dt) {
//   Vec3<Double> mid_x = mid_point(x, x0);
//   Vec3<Double> mid_u(0.5 * (u[0] + u0[0]), 0.5 * (u[1] + u0[1]),
//                      0.5 * (u[2] + u0[2]));
//   // Double gamma = Gamma(mid_x, mid_u, mid_cell, grid);
//   // std::cout << "Gamma is " << gamma.x() << std::endl;
//   // Double u_0 = 0.5 * u0_energy(x0, u0, metric) + 0.5 * u0_energy(x, u,
//   // metric);
//   // std::cout << "u0 is " << u_0.x() << std::endl;

//   if (n == 0) {
//     Double result = 0.0;
//     // Double gamma = Gamma(mid_x, mid_u, metric);
//     result +=
//         0.5 * D<1>(metric.g00)(mid_x[0], mid_x[1], mid_x[2]) * cache.gamma;
//     result += 0.5 * D<1>(metric.g11)(mid_x[0], mid_x[1], mid_x[2]) * mid_u[0] *
//               mid_u[0] / (cache.gamma * cache.g11 * cache.g11);
//     result += 0.5 * D<1>(metric.g22)(mid_x[0], mid_x[1], mid_x[2]) * mid_u[1] *
//               mid_u[1] / (cache.gamma * cache.g22 * cache.g22);
//     result += 0.5 * D<1>(metric.g33)(mid_x[0], mid_x[1], mid_x[2]) * mid_u[2] *
//               mid_u[2] / (cache.gamma * cache.g33 * cache.g33);
//     return u[0] - u0[0] - result * dt;
//   } else if (n == 1) {
//     // Double gamma = Gamma(mid_x, mid_u, metric);
//     return u[1] - u0[1] -
//            0.5 * dt * D<2>(metric.g33)(mid_x[0], mid_x[1], mid_x[2]) *
//                mid_u[2] * mid_u[2] / (cache.gamma * cache.g33 * cache.g33);
//   } else {
//     return 0.0;
//   }
// }

// template <typename Double>
// Double
// Fx(int n, const Vec3<double>& u0, const Vec3<double>& x0, const Vec3<Double>& u,
//    const Vec3<Double>& x, const MetricCache<Double>& cache, double dt) {
//   // Vec3<Double> mid_x = mid_point(x, x0);
//   // Vec3<Double> mid_u(0.5 * (u[0] + u0[0]), 0.5 * (u[1] + u0[1]),
//   // 0.5 * (u[2] + u0[2]));

//   // Double gamma = Gamma(mid_x, mid_u, metric);
//   if (n == 0) {
//     return x[0] - x0[0] - dt * u[0] / (cache.gamma * cache.g11);
//   } else if (n == 1) {
//     return x[1] - x0[1] - dt * u[1] / (cache.gamma * cache.g22);
//   } else if (n == 2) {
//     return x[2] - x0[2] - dt * u[2] / (cache.gamma * cache.g33);
//   } else {
//     return 0.0;
//   }
// }

// template <typename Metric, typename Double>
// int
// iterate_newton(Particle& p, const Metric& metric, MetricCache<Double>& cache,
//                double dt) {
//   const int max_steps = 100;
//   const double tolerance = 1.0e-10;

//   Particle p0 = p;

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
//     Vec3<Double> mid_x = mid_point(x, p0.x);
//     Vec3<Double> mid_u = mid_point(u, p0.u);
//     cache.g11 = metric.g11(mid_x[0], mid_x[1], mid_x[2]);
//     cache.g22 = metric.g22(mid_x[0], mid_x[1], mid_x[2]);
//     cache.g33 = metric.g33(mid_x[0], mid_x[1], mid_x[2]);
//     cache.gamma = Gamma(mid_x, mid_u, metric, cache);
//     // Compute the residual function
//     f[0] = Fx(0, p0.u, p0.x, u, x, cache, dt);
//     f[1] = Fx(1, p0.u, p0.x, u, x, cache, dt);
//     f[2] = Fx(2, p0.u, p0.x, u, x, cache, dt);
//     f[3] = Fu(0, p0.u, p0.x, u, x, metric, cache, dt);
//     f[4] = Fu(1, p0.u, p0.x, u, x, metric, cache, dt);
//     f[5] = Fu(2, p0.u, p0.x, u, x, metric, cache, dt);
//     // Fill the Jacobi matrix
//     std::array<std::array<double, 6>, 6> J;
//     std::array<double, 6> F, new_f;
//     for (int m = 0; m < 6; m++) {
//       F[m] = f[m].x();
//       for (int n = 0; n < 6; n++) {
//         J[n][m] = f[n].d(m);
//       }
//     }
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
  ptc.x[1] = 3.1415926535898 * 0.5;
  double L = sqrt(4.0);
  // double L = 100.0;
  ptc.u[2] = L;
  // ptc.x[0] = L * L + sqrt(L * L - 3.0) * L;
  ptc.x[0] = 1.5 + 0.00000000001;
  ptc.is_photon = true;

  Schwarzschild<var> metric(0.5);

  const double dt = 0.1;
  // Particle ptc_initial = ptc;
  std::ofstream out("schwarzschild-isco-0.1.txt", std::ofstream::out);
  out << std::setprecision(std::numeric_limits<double>::digits10 + 2);
  timer::stamp();

  // std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 2);
  var u_0 = u0_energy(ptc.x, ptc.u, metric, ptc.is_photon);
  // auto pos = grid.mesh().pos_particle(p.cell, p.x);
  std::cout << "Initial: " << ptc << ", u_0 = " << u_0.x() << std::endl;
  const int N = 1000;
  for (int n = 0; n < N; n++) {
    // std::cout << "At timestep " << n << std::endl;

    // double u_0 = u0_energy(ptc.x, ptc.u, metric_num);
    // auto pos = grid.mesh().pos_particle(p.cell, p.x);
    // auto pos = ptc.x;
    // std::cout << pos << " " << pos.x * sin(pos.y) << " " << ptc.u << " " << u_0 << std::endl;

    // if (n % 10 == 0)
    // out << pos.x << ", " << pos.y << ", " << pos.z << ", " << u_0 <<
    // std::endl;

    // if (iterate_newton(ptc, metric, dt, type) == 100) {
    //   std::cout << "Iteration reached end without converging!" << std::endl;
    //   break;
    // }
    // iterate_rk4(ptc, metric, dt);
    iterate_newton(ptc, metric, dt, type);
    std::cout << ptc << std::endl;
  }
  var u_1 = u0_energy(ptc.x, ptc.u, metric, ptc.is_photon);
  // auto pos = grid.mesh().pos_particle(p.cell, p.x);
  std::cout << "Final: " << ptc << ", u_0 = " << u_1.x()
            << ", \\Delta u_0/u_0 = " << (u_1.x() - u_0.x()) / u_0.x() << std::endl;
  std::cout << "Evolved for " << N << " timesteps." << std::endl;
  timer::show_duration_since_stamp("Evolution time", "ms");
  out.close();

  return 0;
}
