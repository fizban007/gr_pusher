#include "CudaLE.h"
#include "fadiff.h"
#include "utils/timer.h"
#include "vec3.h"
#include "ptc.h"
// #include "metrics/metric_cks.h"
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


// typedef fadbad::F<double, 6> var;

int
main(int argc, char* argv[]) {
  // SolverType type = SolverType::implicit;
  SolverType type = SolverType::hamiltonian;
  Particle<var> ptc;
  double m = 1.0;
  double a = 1.0;
  Cartesian_KS<var> metric(m, a);
  // ptc.x[1] = 3.1415926535898 * 0.5;
  // double L = sqrt(4.0) + 1.0e-12;
  // ptc.x[1] = L * L + sqrt(L * L - 3.0) * L;

  // L = 0
  // ptc.x[0] = 2.414213562373095;
  // ptc.x[1] = -1.0;
  // ptc.u[0] = 0.6306019374818707;
  // ptc.u[1] = -0.26120387496374076;
  // ptc.u[2] = 0.511081084529394;

  // L = 1.36
  ptc.x[0] = 1.8;
  ptc.x[1] = -1.0;
  ptc.u[0] = 0.6902515723270438;
  ptc.u[1] = -0.2575471698113204;
  ptc.u[2] = 0.33166247903554;
  ptc.is_photon = true;

  // Cartesian_KS<var> metric(1.0, 1.0);
  // Cartesian_KS<double> metric_num(1.0, 1.0);
  const double dt = 0.001;
  // Particle ptc_initial = ptc;
  std::ofstream out("cks-L1.36-t-3-ham.txt", std::ofstream::out);
  out << std::setprecision(std::numeric_limits<double>::digits10 + 2);
  timer::stamp();

  var u_0 = u0_energy(ptc.x, ptc.u, metric, ptc.is_photon);
  double e_initial = -u_0.x();
  std::cout << "Initial: " << ptc << ", u_0 = " << u_0.x() << std::endl;

  int N = 100000;
  for (int n = 0; n < N; n++) {
    std::cout << "At timestep " << n << std::endl;

    // auto pos = grid.mesh().pos_particle(p.cell, p.x);
    // auto pos = ptc.x;
    // std::cout << ptc << " " << u_0.x()
              // << std::endl;

    int iter;
    if ((iter = iterate_newton(ptc, metric, dt, type)) == -1) {
      std::cout << "Iteration reached end without converging!" << std::endl;
      break;
    }
    u_0 = u0_energy(ptc.x, ptc.u, metric, ptc.is_photon);
    var uu = metric.inv_g00(ptc.x[0], ptc.x[1], ptc.x[2]) * u_0 * u_0 +
             metric.inv_g11(ptc.x[0], ptc.x[1], ptc.x[2]) * ptc.u[0] * ptc.u[0] +
             metric.inv_g22(ptc.x[0], ptc.x[1], ptc.x[2]) * ptc.u[1] * ptc.u[1] +
             metric.inv_g33(ptc.x[0], ptc.x[1], ptc.x[2]) * ptc.u[2] * ptc.u[2] +
             2.0 * metric.inv_g01(ptc.x[0], ptc.x[1], ptc.x[2]) * u_0 * ptc.u[0] +
             2.0 * metric.inv_g02(ptc.x[0], ptc.x[1], ptc.x[2]) * u_0 * ptc.u[1] +
             2.0 * metric.inv_g03(ptc.x[0], ptc.x[1], ptc.x[2]) * u_0 * ptc.u[2] +
             2.0 * metric.inv_g12(ptc.x[0], ptc.x[1], ptc.x[2]) * ptc.u[0] * ptc.u[1] +
             2.0 * metric.inv_g13(ptc.x[0], ptc.x[1], ptc.x[2]) * ptc.u[0] * ptc.u[2] +
             2.0 * metric.inv_g23(ptc.x[0], ptc.x[1], ptc.x[2]) * ptc.u[1] * ptc.u[2];
    double r = metric.r(ptc.x[0], ptc.x[1], ptc.x[2]).x();
    std::cout << "E diff = " << std::abs(u_0.x() + e_initial) / e_initial << ", u2 = " << uu.x() << ", r = " << std::abs(r - 1.8) << ", in " << iter << " iterations" << std::endl;
    if (n % 1000 == 0)
      out << std::abs(u_0.x() + e_initial) / e_initial << ", " << uu.x() << ", " << std::abs(r - 1.8) << ", " << iter << std::endl;

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
