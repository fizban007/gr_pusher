#include <iostream>
#include "fadiff.h"
#include "algorithms/interp_template.h"
#include "catch.hpp"
#include "grid.h"

using namespace Aperture;
using namespace fadbad;

template <typename Double>
Double func(Double x, Double y) {
  Double result = 0.0;
  for (int i = 2; i < 5; i++) {
    for (int j = 1; j < 4; j++) {
      result += interp_cell(x, i - 4, 1.0, 1) * interp_cell(y, j - 3, 1.0, 1) * (double)i * (double)j;
    }
  }
  return result;
}

TEST_CASE("Derivative with respect to interpolation", "[derivative]") {
  F<double> x, f;
  x = 0.3;
  x.diff(0, 1);
  f = interp_cell(x, 1, 1.0);
  // f.diff(0, 1);
  double fval = f.x();
  double dfdx = f.d(0);

  std::cout << "f(x) = " << fval << std::endl;
  std::cout << "dfdx(x) = " << dfdx << std::endl;
}

TEST_CASE("More complex derivative with respect to interpolation", "[derivative]") {
  F<double> x, y, f;
  x = 0.3;
  x.diff(0, 2);
  y = 0.4;
  y.diff(1, 2);
  f = func(x, y);

  double fval = f.x();
  double dfdx = f.d(0);
  double dfdy = f.d(1);

  std::cout << "f(x) = " << fval << std::endl;
  std::cout << "dfdx(x) = " << dfdx << std::endl;
  std::cout << "dfdy(x) = " << dfdy << std::endl;
}

TEST_CASE("Grid interpolation of alpha", "[grid][derivative]") {
  F<double> x1, x2, x3, f;

  Grid g({"DIM1 64 1.0 5.00 2",
          "DIM2 64 0.0 3.14 2",
          "DIM3 1  0.0 0.01 0"});
  double a = 0.99, rg = 2.4;
  auto m = metric::metric_kerr_schild(rg, a);
  g.setup_metric(m, g);

  x1 = 0.3 * g.mesh().delta[0]; x1.diff(0, 3);
  x2 = 0.4 * g.mesh().delta[1]; x2.diff(1, 3);
  x3 = 0.3 * g.mesh().delta[2]; x3.diff(2, 3);

  f = g.beta(0, 5 + 32 * 64, x1, x2, x3);
  double fval = f.x();
  double dfdx1 = f.d(0);
  double dfdx2 = f.d(1);
  double dfdx3 = f.d(2);

  std::cout << "f(x) = " << fval << std::endl;
  std::cout << "dfdx1(x) = " << dfdx1 << std::endl;
  std::cout << "dfdx2(x) = " << dfdx2 << std::endl;
  std::cout << "dfdx3(x) = " << dfdx3 << std::endl;
}
