#include "quadmesh.h"
#include "catch.hpp"

using namespace Aperture;

TEST_CASE("Default constructor", "[quadmesh]") {
  Quadmesh mesh;
  for (int i = 0; i < 3; i++) {
    REQUIRE(mesh.dims[i] == 1);
    REQUIRE(mesh.guard[i] == 0);
    REQUIRE(mesh.delta[i] == Approx(1.0));
    REQUIRE(mesh.lower[i] == Approx(0.0));
    REQUIRE(mesh.sizes[i] == Approx(0.0));
  }
  REQUIRE(mesh.dimension == 1);

  Quadmesh mesh2 (10, 20, -1);
  REQUIRE(mesh2.dims[0] == 10);
  REQUIRE(mesh2.dims[1] == 20);
  REQUIRE(mesh2.dims[2] == 1);
  REQUIRE(mesh2.dimension == 2);
}

TEST_CASE("Testing pos_3d", "[quadmesh]") {
  Quadmesh mesh (20, 20);
  mesh.delta[0] = mesh.delta[1] = 0.1;
  auto v = mesh.pos_3d(210, Stagger_t("001"));
  REQUIRE(v.x == Approx(1.1));
  REQUIRE(v.y == Approx(1.05));
}

TEST_CASE("Dimension", "[quadmesh]") {
  Quadmesh m0 (20, 20, 20);
  REQUIRE(m0.dim() == 3);

  Quadmesh m1 (20, 20);
  REQUIRE(m1.dim() == 2);

  Quadmesh m2 (20);
  REQUIRE(m2.dim() == 1);
}
