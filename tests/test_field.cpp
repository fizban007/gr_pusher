#include <iostream>
#include "grid.h"
#include "fields.h"
#include "metrics.h"
#include "CudaLE.h"
#include "catch.hpp"

using namespace std;
using namespace CudaLE;
using namespace CudaLE::placeholders::spherical;
using namespace Aperture;

Approx zero = Approx(0.0);

TEST_CASE("Initialize", "[scalarfield]") {
  Grid g(10, 10);
  ScalarField<Scalar> field(g);
  field.initialize();
  // This is equivalent to the following:
  // field.assign(0.0);

  REQUIRE(field(0, 0) == zero);
  REQUIRE(field(9, 9) == zero);
  REQUIRE(field.data().size() == 100);
  REQUIRE(field.grid_size() == 100);
}

TEST_CASE("Nontrivial initialize", "[scalarfield]") {
  Grid g({"DIM1 32 1.0 5.00 2",
          "DIM2 32 0.0 3.14 2",
          "DIM3 32 0.0 6.28 2"});
  ScalarField<Scalar> f(g);
  f.initialize( sin(_theta) * cos(_phi) / (_r * _r) );

  REQUIRE(f(10, 10, 10) - std::sin(8.5 * 3.14 / 32.0) * std::cos(8.5 * 6.28 / 32.0) /
                  ((1.0 + 8.5 * 5.0 / 32.0) * (1.0 + 8.5 * 5.0 / 32.0)) == zero);
}

TEST_CASE("Copy", "[scalarfield]") {
  Grid g(10, 10);
  ScalarField<Scalar> f1(g);
  ScalarField<Scalar> f2(g);
  f1.assign(3.0);
  f2.copyFrom(f1);

  REQUIRE(f2(0, 0) == Approx(3.0));
  REQUIRE(f2(9, 9) == Approx(3.0));
  REQUIRE(f2.data().size() == 100);
}

TEST_CASE("Multiply", "[scalarfield]") {
  Grid g(10, 10);
  ScalarField<Scalar> f1(g);
  ScalarField<Scalar> f2(g);
  f1.assign(3.0);
  f2.assign(2.0);

  f1.multiplyBy(3.0);

  REQUIRE(f1(0, 0) - 9.0 == zero);
  REQUIRE(f1(9, 9) - 9.0 == zero);

  f1.multiplyBy(f2);

  REQUIRE(f1(0, 0) - 18.0 == zero);
  REQUIRE(f1(9, 9) - 18.0 == zero);
}

TEST_CASE("Plus", "[scalarfield]") {
  Grid g(10, 10);
  ScalarField<Scalar> f1(g);
  ScalarField<Scalar> f2(g);
  f1.assign(3.0);
  f2.assign(2.0);

  f1.addBy(3.0);

  REQUIRE(f1(0, 0) - 6.0 == zero);
  REQUIRE(f1(9, 9) - 6.0 == zero);

  f1.addBy(f2);

  REQUIRE(f1(0, 0) - 8.0 == zero);
  REQUIRE(f1(9, 9) - 8.0 == zero);
}

TEST_CASE("Chaining operations", "[scalarfield]") {
  Grid g(10, 10);
  ScalarField<Scalar> f1(g);
  ScalarField<Scalar> f2(g);
  f1.assign(3.0);
  f2.assign(2.0);

  // This statement will first add f1 by 3.0 and then multiply it by
  // 1.2, finally add to it the values of f2
  f1.addBy(3.0).multiplyBy(1.2).addBy(f2);

  REQUIRE(f1(0, 0) - 9.2 == zero);
  REQUIRE(f1(9, 9) - 9.2 == zero);
}

TEST_CASE("Initialize Vector field", "[vectorfield]") {
  Grid g({"DIM1 32 1.0 5.00 2",
          "DIM2 32 0.0 3.14 2",
          "DIM3 32 0.0 6.28 2"});
  VectorField<Scalar> f(g);
  f.set_field_type(FieldType::E);
  f.initialize(1, sin(_theta) * cos(_phi) / (_r * _r) );

  REQUIRE(f(1, 10, 10, 10) - std::sin(9.0 * 3.14 / 32.0) * std::cos(8.5 * 6.28 / 32.0) /
                  ((1.0 + 8.5 * 5.0 / 32.0) * (1.0 + 8.5 * 5.0 / 32.0)) == zero);
}

TEST_CASE("Interpolate", "[vectorfield]") {
  Grid g({"DIM1 32 1.0 5.00 2",
          "DIM2 32 0.0 3.14 2",
          "DIM3 32 0.0 6.28 2"});
  VectorField<Scalar> f(g);
  f.set_field_type(FieldType::E);
  f.assign(0.0);
  f.initialize(0, _r);
  f.initialize(1, _theta);

  auto v = f.interpolate(Vec3<int>(10, 5, 5), Vec3<Pos_t>(0.3, 0.6, 0.0), 1);
  // ASSERT_FLOAT_EQ(v[2], 0.0);
  REQUIRE(v[0] - g.mesh().pos(0, 10, 0.3) == zero);
  REQUIRE(v[1] - g.mesh().pos(1, 5, 0.6) == zero);
}

TEST_CASE("Recenter", "[vectorfield]") {
  Grid g({"DIM1 32 1.0 5.00 2",
          "DIM2 32 0.0 3.14 2",
          "DIM3 32 0.0 6.28 2"});
  VectorField<Scalar> f(g), f2(g);
  f.set_field_type(FieldType::E);
  f.assign(0.0);
  f.initialize(0, _r);
  f.initialize(1, _theta);

  f.recenter(f2);

  REQUIRE(f2(0, 10, 5, 5) - g.mesh().pos(0, 10, 0.5) == zero);
  REQUIRE(f2(1, 10, 5, 5) - g.mesh().pos(1, 5, 0.5) == zero);
}
