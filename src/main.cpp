#include <iostream>
#include <Eigen/Dense>
#include "fadiff.h"
#include "grid.h"
#include "metrics.h"
#include "vec3.h"

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
  // std::cout << "Mid point is (" << result[0].x() << ", " << result[1].x() << ", " << result[2].x() << ")" << std::endl;
  return result;
}

template <typename Double>
Double
Gamma(const Vec3<Double>& x, const Vec3<Double>& u, int cell,
      const Grid& grid) {
  Double result = 1.0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      result += grid.inv_metric(i, j, cell, x) * u[i] * u[j];
  // std::cout << "alpha is " << grid.alpha(cell, x).x() << std::endl;
  result = sqrt(result) / grid.alpha(cell, x);
  return result;
}

template <typename Double>
Double
Fu(int n, const Vec3<double>& u0, const Vec3<double>& x0, int c0,
   const Vec3<Double>& u, const Vec3<Double>& x, int cell, const Grid& grid,
   double dt) {
  int mid_cell = 0;
  Vec3<Double> mid_x = mid_point(x, x0, cell, c0, mid_cell, grid);
  Vec3<Double> mid_u(0.5 * (u[0] + u0[0]), 0.5 * (u[1] + u0[1]),
                     0.5 * (u[2] + u0[2]));
  Double gamma = Gamma(mid_x, mid_u, mid_cell, grid);
  // std::cout << "Gamma is " << gamma.x() << std::endl;
  Double u_0 = grid.beta(0, mid_cell, mid_x) * mid_u[0] + grid.beta(1, mid_cell, mid_x) * mid_u[1] + grid.beta(2, mid_cell, mid_x) * mid_u[2] - grid.alpha(mid_cell, mid_x) * grid.alpha(mid_cell, mid_x) * gamma;
  // std::cout << "u0 is " << u_0.x() << std::endl;

  Double conn = 0.0;
  for (int i = 1; i < 4; i++) {
    for (int j = 1; j < 4; j++) {
      if (grid.conn_mask(n, i, j))
        conn += grid.connection(n, i, j, mid_cell, mid_x) *
                mid_u[i - 1] * mid_u[j - 1];
    }
  }
  for (int i = 0; i < 3; i++) {
    if (grid.conn_mask(n, i + 1, 0))
      conn += grid.connection(n, i + 1, 0, mid_cell, mid_x) * mid_u[i] * u_0;
    if (grid.conn_mask(n, 0, i + 1))
      conn += grid.connection(n, 0, i + 1, mid_cell, mid_x) * mid_u[i] * u_0;
  }
  conn += grid.connection(n, 0, 0, mid_cell, mid_x) * u_0 * u_0;
  // return u[n] - u0[n] - dt * conn / (Gamma(x, u, cell, grid) + Gamma(x0, u0, c0, grid));
  return u[n] - u0[n] - 0.5 * dt * conn / Gamma(mid_x, mid_u, mid_cell, grid);
}

template <typename Double>
Double
Fx(int n, const Vec3<double>& u0, const Vec3<double>& x0, int c0,
   const Vec3<Double>& u, const Vec3<Double>& x, int c, const Grid& grid,
   double dt) {
  Double result = 0.0;
  int mid_cell = 0;
  auto mid_x = mid_point(x, x0, c, c0, mid_cell, grid);
  Vec3<Double> mid_u(0.5 * (u[0] + u0[0]), 0.5 * (u[1] + u0[1]),
                     0.5 * (u[2] + u0[2]));
  auto cell = grid.mesh().get_cell_3d(c);
  auto cell_0 = grid.mesh().get_cell_3d(c0);
  for (int j = 0; j < 3; j++) {
    result -= dt * grid.inv_metric(n, j, mid_cell, mid_x) * (u[j] +  u0[j]);
  }
  // std::cout << "gamma u2 is (" << (grid.inv_metric(n, 1, mid_cell, mid_x) * (u0[1] + u[1])).x() << ")" << std::endl;
  // result /= (Gamma(x, u, c, grid) + Gamma(x0, u0, c0, grid));
  result /= 2.0 * Gamma(mid_x, mid_u, mid_cell, grid);
  result += x[n] - x0[n] - (cell_0[n] - cell[n]) * grid.mesh().delta[n] + dt * grid.beta(n, mid_cell, mid_x);
  return result;
}

int iterate_newton(Particle& p, const Grid& grid, double dt) {
  const int max_steps = 100;
  const double tolerance = 1.0e-10;
  typedef F<double, 6> Var;
  Particle p0 = p;
  // std::cout << "Initially, " << p0.x << " " << p0.u << " " << p0.cell << std::endl;
  for (int i = 0; i < max_steps; i++) {
    // std::cout << "at step " << i << std::endl;
    // Initialize guess solution
    Vec3<Var> x(p.x[0], p.x[1], p.x[2]);
    Vec3<Var> u(p.u[0], p.u[1], p.u[2]);
    x[0].diff(0);
    x[1].diff(1);
    x[2].diff(2);
    u[0].diff(3);
    u[1].diff(4);
    u[2].diff(5);
    std::array<Var, 6> f;
    int c = p.cell;
    // Compute the residual function
    f[0] = Fx(0, p0.u, p0.x, p0.cell, u, x, c, grid, dt);
    f[1] = Fx(1, p0.u, p0.x, p0.cell, u, x, c, grid, dt);
    f[2] = Fx(2, p0.u, p0.x, p0.cell, u, x, c, grid, dt);
    f[3] = Fu(0, p0.u, p0.x, p0.cell, u, x, c, grid, dt);
    f[4] = Fu(1, p0.u, p0.x, p0.cell, u, x, c, grid, dt);
    f[5] = Fu(2, p0.u, p0.x, p0.cell, u, x, c, grid, dt);
    // Fill the Jacobi matrix
    Eigen::MatrixXd J(6, 6);
    Eigen::VectorXd F(6);
    for (int m = 0; m < 6; m++) {
      for (int n = 0; n < 6; n++) {
        J(n, m) = f[n].d(m);
      }
      F(m) = f[m].x();
    }
    // std::cout << J << std::endl;
    // std::cout << F << std::endl;
    // auto J_inv = J.inverse();
    Eigen::VectorXd new_f = J.colPivHouseholderQr().solve(F);
    // std::cout << new_f << std::endl;
    p.x[0] -= new_f[0];
    p.x[1] -= new_f[1];
    p.x[2] -= new_f[2];
    p.u[0] -= new_f[3];
    p.u[1] -= new_f[4];
    p.u[2] -= new_f[5];
    auto new_cell = grid.mesh().get_cell_3d(p.cell);
    for (int j = 0; j < grid.dim(); j++) {
      if (p.x[j] >= grid.mesh().delta[j]) {
        new_cell[j] += 1;
        p.x[j] -= grid.mesh().delta[j];
      } else if (p.x[j] < 0.0) {
        new_cell[j] -= 1;
        p.x[j] += grid.mesh().delta[j];
      }
    }
    p.cell = grid.mesh().get_idx(new_cell[0], new_cell[1], new_cell[2]);
    // if (sqrt(new_f.dot(new_f)) <= tolerance) break;
    // std::cout << "Norm is " << new_f.norm() << std::endl;
    // double u_0 = grid.beta(0, p.cell, p.x) * p.u[0] + grid.beta(1, p.cell, p.x) * p.u[1] + grid.beta(2, p.cell, p.x) * p.u[2] - grid.alpha(p.cell, p.x) * grid.alpha(p.cell, p.x) * Gamma(p.x, p.u, p.cell, grid);
    // std::cout << p.x << " " << p.u << " " << grid.mesh().get_cell_3d(p.cell) << " " << u_0 << std::endl;
    if (new_f.norm() <= tolerance) break;
    if (i == max_steps - 1) return 1;
  }
  double u_0 = grid.beta(0, p.cell, p.x) * p.u[0] + grid.beta(1, p.cell, p.x) * p.u[1] + grid.beta(2, p.cell, p.x) * p.u[2] - grid.alpha(p.cell, p.x) * grid.alpha(p.cell, p.x) * Gamma(p.x, p.u, p.cell, grid);
  auto pos = grid.mesh().pos_particle(p.cell, p.x);
  std::cout << pos << " " << pos.x * sin(pos.y) <<  " " << p.u << " " << u_0 << std::endl;
  return 0;
}

int
main(int argc, char* argv[]) {
  ////////////////////////////////////////////////////////////////////////////////
  ///  Setting up the grid and metric we are going to use.
  ////////////////////////////////////////////////////////////////////////////////

  // The convention for grid parsing is as follows:
  // Number of cells, Starting coordinate, Total length, Number of guard cells
  Grid grid(
      {"DIM1 256 1.0 5.00 1", "DIM2 256 0.0 3.14 1", "DIM3 256 0.0 6.28 1"});
      // {"DIM1 256 1.0 5.00 1", "DIM2 256 0.0 3.14 1", "DIM3 1 0.0 6.28 0"});
  auto m = metric::metric_spherical();
  //     {"DIM1 256 0.0 2.00 1", "DIM2 256 0.0 3.14 1", "DIM3 256 0.0 6.28 1"});
  // auto m = metric::metric_log_spherical();

  // This line will cache all the necessary quantities on the grid, and save
  // them into a file. If a file with the same grid/metric configuration is
  // present, it reads from the file instead.
  grid.setup_metric(m, grid);

  ////////////////////////////////////////////////////////////////////////////////
  ///  Initialize particle initial condition.
  ////////////////////////////////////////////////////////////////////////////////
  Particle p;
  // p.cell = 128 + 129 * 258 + 24 * 258 * 258;
  p.cell = 128 + 129 * 258 + 24 * 258 * 258;
  p.u[1] = -1.0;
  // p.x[0] = 0.5 * grid.mesh().delta[0];
  // p.x[1] = 0.5 * grid.mesh().delta[1];
  // p.x[2] = 0.5 * grid.mesh().delta[2];

  ////////////////////////////////////////////////////////////////////////////////
  ///  Set up iteration
  ////////////////////////////////////////////////////////////////////////////////
  const double dt = 0.01;
  for (int n = 0; n < 1000; n++) {
    std::cout << "At timestep " << n << std::endl;
    if (iterate_newton(p, grid, dt) == 1) break;
  }

  return 0;
}
