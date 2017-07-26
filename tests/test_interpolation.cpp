#include "algorithms/interpolation.h"
#include "algorithms/interp_template.h"
#include "fields.h"
#include "catch.hpp"
#include "CudaLE.h"

using namespace Aperture;

TEST_CASE("Approx(0.0)th order", "[interpolation]") {
  Interpolator interp(0);

  REQUIRE(interp.interp_cell(0.6, 5, 6) - 0.0 == Approx(0.0));
  REQUIRE(interp.interp_cell(0.6, 5, 5) - 1.0 == Approx(0.0));
  REQUIRE(interp.interp_cell(0.1, 5, 4) - 0.0 == Approx(0.0));

  REQUIRE(interp.interp_cell(0.3, 5, 4, 1) - 1.0 == Approx(0.0));
}

TEST_CASE("Linear order", "[interpolation]") {
  Interpolator interp(1);

  REQUIRE(interp.interp_cell(0.6, 5, 6) - 0.1 == Approx(0.0));
  REQUIRE(interp.interp_cell(0.6, 5, 5) - 0.9 == Approx(0.0));
  REQUIRE(interp.interp_cell(0.1, 5, 4) - 0.4 == Approx(0.0));

  REQUIRE(interp.interp_cell(0.6, 5, 5, 1) - 0.6 == Approx(0.0));
  REQUIRE(interp.interp_cell(0.6, 5, 4, 1) - 0.4 == Approx(0.0));
}

TEST_CASE("New interpolation", "[interp_template]") {
  REQUIRE(interp_cell(0.6, 1, 1.0) == Approx(0.1));
  REQUIRE(interp_cell(0.6, 0, 1.0) == Approx(0.9));
  REQUIRE(interp_cell(0.1, -1, 1.0) == Approx(0.4));

  REQUIRE(interp_cell(0.6, 0, 1.0, 1) == Approx(0.6));
  REQUIRE(interp_cell(0.6, -1, 1.0, 1) == Approx(0.4));
}

TEST_CASE("Interpolating field", "[interp_template][fields]") {
  Grid grid(
      {"DIM1 256 1.0 5.00 1", "DIM2 256 0.0 3.14 1", "DIM3 1 0.0 6.28 0"});

  auto m = metric::metric_spherical();
  grid.setup_metric(m, grid);

  VectorField<double> E(grid);
  double Ez = 10.0;
  using CudaLE::placeholders::spherical::_theta;
  E.initialize(0, Ez * sin(_theta));
  E.initialize(1, -Ez * cos(_theta));

  std::cout << E.interpolate(0, 127 + 129 * 258, Vec3<double>(0.0, 0.0, 0.0)) << std::endl;

}
