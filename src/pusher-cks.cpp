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
  Particle<var> ptc;
  double m = 1.0;
  double a = 0.0;
  Cartesian_KS<var> metric(m, a);
  // ptc.x[1] = 3.1415926535898 * 0.5;
  // double L = sqrt(4.0) + 1.0e-12;
  ptc.u[0] = 3.0;
  // ptc.x[1] = L * L + sqrt(L * L - 3.0) * L;
  ptc.x[1] = 5.0;

  // Cartesian_KS<var> metric(1.0, 1.0);
  // Cartesian_KS<double> metric_num(1.0, 1.0);
  const double dt = 0.01;
  // Particle ptc_initial = ptc;
  std::ofstream out("cks-isco-0.1.txt", std::ofstream::out);
  out << std::setprecision(std::numeric_limits<double>::digits10 + 2);
  timer::stamp();

  for (int n = 0; n < 10000; n++) {
    std::cout << "At timestep " << n << std::endl;

    var u_0 = u0_energy(ptc.x, ptc.u, metric, ptc.is_photon);
    // auto pos = grid.mesh().pos_particle(p.cell, p.x);
    auto pos = ptc.x;
    std::cout << ptc << " " << u_0.x()
              << std::endl;

    // if (n % 10 == 0)
    //   out << pos.x << ", " << pos.y << ", " << pos.z << ", " << u_0
    //       << std::endl;

    if (iterate_newton(ptc, metric, dt) == 1) {
      std::cout << "Iteration reached end without converging!" << std::endl;
      break;
    }
  }
  timer::show_duration_since_stamp("Evolution time", "ms");
  out.close();

  return 0;
}
