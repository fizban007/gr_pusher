#include "fadiff.h"
#include "grid.h"
#include "metrics.h"
#include "vec3.h"
#include <iostream>

using namespace Aperture;
using namespace fadbad;

struct Particle {
  Vec3<double> x = {0.0, 0.0, 0.0};
  Vec3<double> u = {0.0, 0.0, 0.0};
  int cell = 0;
};  // ----- end of struct Particle -----

template <typename Double>
Vec3<Double>
mid_point(const Vec3<Double>& x, const Vec3<double>& x0, int c, int c0,
          int& c_result, const Grid& grid) {
  // TODO: check the correctness of this function
  Vec3<Double> result{0.0, 0.0, 0.0};
  if (c == c0) {
    c_result = c0;
    for (int i = 0; i < 3; i++) result[i] = 0.5 * (x[i] + x0[i]);
  } else {
    auto cell = grid.mesh().get_cell_3d(c);
    auto cell_0 = grid.mesh().get_cell_3d(c0);
    for (int i = 0; i < 3; i++) {
      result[i] =
          0.5 * (x[i] + x0[i] + (cell_0[i] - cell[i]) * grid.mesh().delta[i]);
      if (result[i] >= grid.mesh().delta[i]) {
        result[i] -= grid.mesh().delta[i];
        cell[i] += 1;
      } else if (result[i] < 0.0) {
        result[i] += grid.mesh().delta[i];
        cell[i] -= 1;
      }
    }
    c_result = grid.mesh().get_idx(cell[0], cell[1], cell[2]);
  }
  return result;
}

template <typename Double>
Double
Gamma(const Vec3<Double>& x, const Vec3<Double>& u, int cell,
      const Grid& grid) {
  Double result = 1.0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      result += grid.inv_metric(i, j, cell, x[0], x[1], x[2]) * u[i] * u[j];
  result = sqrt(result) / grid.alpha(cell, x[0], x[1], x[2]);
  return result;
}

template <typename Double>
Double
Fu(int n, const Vec3<double>& u0, const Vec3<double>& x0, int c0,
   const Vec3<Double>& u, const Vec3<Double>& x, int cell, const Grid& grid,
   double dt) {
  Double result = 0.0;
  int mid_cell = 0;
  auto mid_x = mid_point(x, x0, cell, c0, mid_cell, grid);
  Vec3<Double> mid_u(0.5 * (u[0] + u0[0]), 0.5 * (u[1] + u0[1]),
                     0.5 * (u[2] + u0[2]));
  Double conn = 0.0;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      conn += grid.connection(n, i, j, mid_cell, mid_x[0], mid_x[1], mid_x[2]) *
              mid_u[i] * mid_u[j];
    }
  }
  return u[n] - u0[n] - 0.5 * dt * conn / Gamma(mid_x, mid_u, mid_cell, grid);
}

int
main(int argc, char* argv[]) {
  ////////////////////////////////////////////////////////////////////////////////
  ///  Setting up the grid and metric we are going to use.
  ////////////////////////////////////////////////////////////////////////////////

  // The convention for grid parsing is as follows:
  // Number of cells, Starting coordinate, Total length, Number of guard cells
  Grid grid(
      {"DIM1 256 1.0 5.00 1", "DIM2 256 0.0 3.14 1", "DIM3 1   0.0 0.01 0"});
  auto m = metric::metric_spherical();

  // This line will cache all the necessary quantities on the grid, and save
  // them into a file. If a file with the same grid/metric configuration is
  // present, it reads from the file instead.
  grid.setup_metric(m, grid);

  ////////////////////////////////////////////////////////////////////////////////
  ///  Initialize particle initial condition.
  ////////////////////////////////////////////////////////////////////////////////
  Particle p;
  p.cell = 128 + 129 * 258;
  p.u[1] = -1.0;

  ////////////////////////////////////////////////////////////////////////////////
  ///  Set up iteration
  ////////////////////////////////////////////////////////////////////////////////

  return 0;
}
