#include "CudaLE.h"
#include "fadiff.h"
#include "solve.h"
#include "utils/timer.h"
#include "vec3.h"
// #include "CudaLE.h"
// #include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace Aperture;
using namespace CudaLE;
// using namespace fadbad;

#define CONNECTION(METRIC, I, A, B, X)                                    \
  METRIC.inv_g ## I ## 0(X[0], X[1], X[2]) * (D<B>( METRIC.g ## A ## 0 )(X[0], X[1], X[2]) + D<A>( METRIC.g ## B ## 0)(X[0], X[1], X[2])) + \
  METRIC.inv_g ## I ## 1(X[0], X[1], X[2]) * (D<B>( METRIC.g ## A ## 1 )(X[0], X[1], X[2]) + D<A>( METRIC.g ## B ## 1)(X[0], X[1], X[2]) - D<1>( METRIC.g ## A ## B)(X[0], X[1], X[2])) + \
  METRIC.inv_g ## I ## 2(X[0], X[1], X[2]) * (D<B>( METRIC.g ## A ## 2 )(X[0], X[1], X[2]) + D<A>( METRIC.g ## B ## 2)(X[0], X[1], X[2]) - D<2>( METRIC.g ## A ## B)(X[0], X[1], X[2])) + \
  METRIC.inv_g ## I ## 3(X[0], X[1], X[2]) * (D<B>( METRIC.g ## A ## 3 )(X[0], X[1], X[2]) + D<A>( METRIC.g ## B ## 3)(X[0], X[1], X[2]) - D<3>( METRIC.g ## A ## B)(X[0], X[1], X[2]))


typedef fadbad::F<double, 6> var;

struct Particle {
  Vec3<double> x = {0.0, 0.0, 0.0};
  Vec3<double> u = {0.0, 0.0, 0.0};
  double u0 = 0.0;
  double e_over_m = -1.0;
};  // ----- end of struct Particle -----

template <typename Data>
struct Cartesian_KS {
  CudaLE::Var<1, Data> _x;
  CudaLE::Var<2, Data> _y;
  CudaLE::Var<3, Data> _z;

  double rg_ = 2.0;
  double a_ = 0.0;

  Cartesian_KS() {}
  Cartesian_KS(double rg, double a) : rg_(rg), a_(a) {}

  // DEFINE_FUNCTOR(g00, (rg / _r - 1.0));
  // DEFINE_FUNCTOR(g11, 1.0 / (1.0 - rg / _r));
  // DEFINE_FUNCTOR(g22, _r* _r);
  // DEFINE_FUNCTOR(g33, _r* _r* square(sin(_theta)));

  DEFINE_FUNCTOR(R2, square(_x) + square(_y) + square(_z));
  DEFINE_FUNCTOR(r,
                 sqrt(0.5 * (R2 - a_ * a_ + sqrt(square(R2 - a_ * a_) +
                                                 4.0 * a_ * a_ * square(_z)))));
  DEFINE_FUNCTOR(l1, (r * _x + a_ * _y) / (square(r) + a_ * a_));
  DEFINE_FUNCTOR(l2, (r * _y - a_ * _x) / (square(r) + a_ * a_));
  DEFINE_FUNCTOR(l3, _z / r);
  DEFINE_FUNCTOR(f, rg_* pow<3>(r) / (pow<4>(r) + a_ * a_ * square(_z)));

  DEFINE_FUNCTOR(g00, -1.0 + f);
  DEFINE_FUNCTOR(g01, f* l1);
  DEFINE_FUNCTOR(g02, f* l2);
  DEFINE_FUNCTOR(g03, f* l3);
  DEFINE_FUNCTOR(g10, f* l1);
  DEFINE_FUNCTOR(g11, 1.0 + f * square(l1));
  DEFINE_FUNCTOR(g12, f* l1* l2);
  DEFINE_FUNCTOR(g13, f* l1* l3);
  DEFINE_FUNCTOR(g20, f* l2);
  DEFINE_FUNCTOR(g21, f* l1* l2);
  DEFINE_FUNCTOR(g22, 1.0 + f * square(l2));
  DEFINE_FUNCTOR(g23, f* l2* l3);
  DEFINE_FUNCTOR(g30, f* l3);
  DEFINE_FUNCTOR(g31, f* l1* l3);
  DEFINE_FUNCTOR(g32, f* l2* l3);
  DEFINE_FUNCTOR(g33, 1.0 + f * square(l3));

  DEFINE_FUNCTOR(inv_g00, -1.0 - f);
  DEFINE_FUNCTOR(inv_g01, f* l1);
  DEFINE_FUNCTOR(inv_g02, f* l2);
  DEFINE_FUNCTOR(inv_g03, f* l3);
  DEFINE_FUNCTOR(inv_g10, f* l1);
  DEFINE_FUNCTOR(inv_g11, 1.0 - f * square(l1));
  DEFINE_FUNCTOR(inv_g12, -1.0 * f * l1 * l2);
  DEFINE_FUNCTOR(inv_g13, -1.0 * f * l1 * l3);
  DEFINE_FUNCTOR(inv_g20, f* l2);
  DEFINE_FUNCTOR(inv_g21, -1.0 * f * l1 * l2);
  DEFINE_FUNCTOR(inv_g22, 1.0 - f * square(l2));
  DEFINE_FUNCTOR(inv_g23, -1.0 * f * l2 * l3);
  DEFINE_FUNCTOR(inv_g30, f* l3);
  DEFINE_FUNCTOR(inv_g31, -1.0 * f * l1 * l3);
  DEFINE_FUNCTOR(inv_g32, -1.0 * f * l2 * l3);
  DEFINE_FUNCTOR(inv_g33, 1.0 - f * square(l3));
};

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
quadratic_solve(const Data& a, const Data& b, const Data& c) {
  return (-b - sqrt(b * b - 4.0 * a * c)) / (2.0 * a);
}

template <typename Data>
Data
Gamma(const Vec3<Data>& x, const Vec3<Data>& u,
      const Cartesian_KS<Data>& metric) {
  // Data result = 1.0;
  // TODO: Make this more general
  // result += u[0] * u[0] / metric.g11(x[0], x[1], x[2]);
  // result += u[1] * u[1] / metric.g22(x[0], x[1], x[2]);
  // result += u[2] * u[2] / metric.g33(x[0], x[1], x[2]);
  // result = sqrt(result) / metric.alpha(x[0], x[1], x[2]);
  Data a = metric.g00(x[0], x[1], x[2]);
  Data b = 2.0 * (metric.g01(x[0], x[1], x[2]) * u[0] + metric.g02(x[0], x[1], x[2]) * u[1] +
                  metric.g03(x[0], x[1], x[2]) * u[2]);
  Data c = 1.0 + metric.g11(x[0], x[1], x[2]) * u[0] * u[0] +
           metric.g22(x[0], x[1], x[2]) * u[1] * u[1] +
           metric.g33(x[0], x[1], x[2]) * u[2] * u[2] +
           2.0 * metric.g12(x[0], x[1], x[2]) * u[0] * u[1] +
           2.0 * metric.g13(x[0], x[1], x[2]) * u[0] * u[2] +
           2.0 * metric.g23(x[0], x[1], x[2]) * u[1] * u[2];
  return quadratic_solve(a, b, c);
}

template <typename Data>
Data
u0_energy(const Vec3<Data>& x, const Vec3<Data>& u,
          const Cartesian_KS<Data>& metric) {
  // return -metric.alpha(x[0], x[1], x[2]) * metric.alpha(x[0], x[1], x[2]) *
  //        Gamma(x, u, metric);
  Data result = Gamma(x, u, metric);
  result = metric.g00(x[0], x[1], x[2]) * result;
  result += metric.g01(x[0], x[1], x[2]) * u[0] + metric.g02(x[0], x[1], x[2]) * u[1] +
            metric.g03(x[0], x[1], x[2]) * u[2];
  return result;
}

// template <typename Data>
// Data
// connection(int i, int a, int b, const Vec3<Data>& x,
//            const Cartesian_KS<Data>& metric) {
//   Data result;
//   if (i == 1) {

//   }
//   return result;
// }

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
  // std::cout << "u0 is " << u_0.x() << std::endl;
  Data result = 0.0;
  Data gamma = Gamma(mid_x, mid_u, metric);

  if (n == 0) {
    result -= 0.5 * CONNECTION(metric, 1, 0, 0, mid_x) * gamma +
              CONNECTION(metric, 1, 0, 1, mid_x) * mid_u[0] +
              CONNECTION(metric, 1, 0, 2, mid_x) * mid_u[1] +
              CONNECTION(metric, 1, 0, 3, mid_x) * mid_u[2] +
              0.5 * CONNECTION(metric, 1, 1, 1, mid_x) * mid_u[0] * mid_u[0] / gamma +
              0.5 * CONNECTION(metric, 1, 2, 2, mid_x) * mid_u[1] * mid_u[1] / gamma +
              0.5 * CONNECTION(metric, 1, 3, 3, mid_x) * mid_u[2] * mid_u[2] / gamma +
              CONNECTION(metric, 1, 1, 2, mid_x) * mid_u[0] * mid_u[1] / gamma +
              CONNECTION(metric, 1, 1, 3, mid_x) * mid_u[0] * mid_u[2] / gamma +
              CONNECTION(metric, 1, 2, 3, mid_x) * mid_u[1] * mid_u[2] / gamma;
  } else if (n == 1) {
    result -= 0.5 * CONNECTION(metric, 2, 0, 0, mid_x) * gamma +
              CONNECTION(metric, 2, 0, 1, mid_x) * mid_u[0] +
              CONNECTION(metric, 2, 0, 2, mid_x) * mid_u[1] +
              CONNECTION(metric, 2, 0, 3, mid_x) * mid_u[2] +
              0.5 * CONNECTION(metric, 2, 1, 1, mid_x) * mid_u[0] * mid_u[0] / gamma +
              0.5 * CONNECTION(metric, 2, 2, 2, mid_x) * mid_u[1] * mid_u[1] / gamma +
              0.5 * CONNECTION(metric, 2, 3, 3, mid_x) * mid_u[2] * mid_u[2] / gamma +
              CONNECTION(metric, 2, 1, 2, mid_x) * mid_u[0] * mid_u[1] / gamma +
              CONNECTION(metric, 2, 1, 3, mid_x) * mid_u[0] * mid_u[2] / gamma +
              CONNECTION(metric, 2, 2, 3, mid_x) * mid_u[1] * mid_u[2] / gamma;
  } else {
    result -= 0.5 * CONNECTION(metric, 3, 0, 0, mid_x) * gamma +
              CONNECTION(metric, 3, 0, 1, mid_x) * mid_u[0] +
              CONNECTION(metric, 3, 0, 2, mid_x) * mid_u[1] +
              CONNECTION(metric, 3, 0, 3, mid_x) * mid_u[2] +
              0.5 * CONNECTION(metric, 3, 1, 1, mid_x) * mid_u[0] * mid_u[0] / gamma +
              0.5 * CONNECTION(metric, 3, 2, 2, mid_x) * mid_u[1] * mid_u[1] / gamma +
              0.5 * CONNECTION(metric, 3, 3, 3, mid_x) * mid_u[2] * mid_u[2] / gamma +
              CONNECTION(metric, 3, 1, 2, mid_x) * mid_u[0] * mid_u[1] / gamma +
              CONNECTION(metric, 3, 1, 3, mid_x) * mid_u[0] * mid_u[2] / gamma +
              CONNECTION(metric, 3, 2, 3, mid_x) * mid_u[1] * mid_u[2] / gamma;
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
           dt * u[0] / gamma;
  } else if (n == 1) {
    return x[1] - x0[1] -
           dt * u[1] / gamma;
  } else if (n == 2) {
    return x[2] - x0[2] -
           dt * u[2] / gamma;
  } else {
    return 0.0;
  }
}

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
    std::cout << "(" << F[0] << ", " << F[1] << ", " << F[2] << ", " << F[3]
              << ", " << F[4] << ", " << F[5] << ")" << std::endl;
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

int
main(int argc, char* argv[]) {
  Particle ptc;
  // ptc.x[1] = 3.1415926535898 * 0.5;
  double L = sqrt(4.0) + 1.0e-12;
  ptc.u[0] = -L;
  ptc.x[1] = L * L + sqrt(L * L - 3.0) * L;

  Cartesian_KS<var> metric(1.0, 0.0);
  Cartesian_KS<double> metric_num(1.0, 0.0);
  const double dt = 0.1;
  // Particle ptc_initial = ptc;
  std::ofstream out("cks-isco-0.1.txt", std::ofstream::out);
  out << std::setprecision(std::numeric_limits<double>::digits10 + 2);
  timer::stamp();

  for (int n = 0; n < 100000; n++) {
    std::cout << "At timestep " << n << std::endl;

    double u_0 = u0_energy(ptc.x, ptc.u, metric_num);
    // auto pos = grid.mesh().pos_particle(p.cell, p.x);
    auto pos = ptc.x;
    std::cout << pos << " " << ptc.u << " " << u_0
              << std::endl;

    if (n % 10 == 0)
      out << pos.x << ", " << pos.y << ", " << pos.z << ", " << u_0
          << std::endl;

    if (iterate_newton(ptc, metric, dt) == 1) {
      std::cout << "Iteration reached end without converging!" << std::endl;
      break;
    }
  }
  timer::show_duration_since_stamp("Evolution time", "ms");
  out.close();

  return 0;
}
