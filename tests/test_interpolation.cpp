#include "algorithms/interpolation.h"
#include "algorithms/interp_template.h"
#include "catch.hpp"

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
