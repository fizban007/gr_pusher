#include "catch.hpp"
#include "CudaLE.h"

// using namespace Aperture;
using namespace CudaLE;
using namespace CudaLE::placeholders;

TEST_CASE("Testing Derivatives", "[CudaLE]") {
  auto f = ConstOp(3.0);
  REQUIRE(D<1>(f) == ZeroOp());

  auto g = placeholders::_2;
  REQUIRE(D<1>(g) == ZeroOp());
  REQUIRE(D<2>(g) == ConstOp(1.0));

  REQUIRE(D<1>(4.0) == ZeroOp());

  REQUIRE(D<1>(ZeroOp()) == ZeroOp());
}

TEST_CASE("Testing Simplify", "[CudaLE]") {
  auto h = ZeroOp{} + ZeroOp{};
  REQUIRE(simplify(h) == ZeroOp{});

  auto g = square(exp(_1) * sin(_2));
  auto f = D<3>(g);
  f.print(); printf("\n");
  D<3>(exp(_1) * sin(_2)).print(); printf("\n");
  // simplify(f).print(); printf("\n");
  // EXPECT_EQ(simplify(f), ZeroOp{});
}
