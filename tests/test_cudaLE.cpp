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
  REQUIRE(D<2>(g) == OneOp());

  REQUIRE(D<1>(4.0) == ZeroOp());

  REQUIRE(D<1>(ZeroOp()) == ZeroOp());
}

TEST_CASE("Testing Simplify", "[CudaLE]") {
  auto h = ZeroOp{} + ZeroOp{};
  REQUIRE(simplify(h) == ZeroOp{});

  auto g = square(exp(_1) * sin(_2));
  auto f = D<3>(g);
  simplify(f).print(); printf("\n");
  simplify(D<3>(exp(_1) * sin(_2))).print(); printf("\n");
  // simplify(f).print(); printf("\n");
  // EXPECT_EQ(simplify(f), ZeroOp{});
}

TEST_CASE("Testing Unary Simplify", "[CudaLE]") {
  auto f = square((_1 * _1) + ZeroOp{});
  simplify(f).print(); printf("\n");
}

TEST_CASE("Testing Derivative Simplify", "[CudaLE]") {
  auto f = D<1>(square(sin(_1) * _2 + 3.0));
  simplify(f).print(); printf("\n");

  auto g = (square(_1 * _1 + 2.0 * 3.0) -
            square(3.0 * sin(_2)) * (_1 * _1 + 4.0)) *
           square(sin(_2)) / (_1 * _1 + 2.0);
  double m_a = 0.6;
  double m_rg = 2.0;
  DEFINE_FUNCTOR(rho2, _1 * _1 + m_a * m_a * square(cos(_2)));
  DEFINE_FUNCTOR(z, _1 * m_rg / rho2);
  DEFINE_FUNCTOR(delta, _1 * _1 + m_a * m_a - m_rg * _1);
  DEFINE_FUNCTOR(g33, (square(_1 * _1 + m_a * m_a) -
                       square(m_a * sin(_2)) * delta) *
                 square(sin(_2)) / rho2);

  // simplify(D<1>(g33)).print(); printf("\n");
  // simplify(D<3>(g33)).print(); printf("\n");
  // simplify(ZeroOp{} / _2).print(); printf("\n");
  std::cout << "Simplify test" << std::endl;
  simplify((ZeroOp{} - ZeroOp{}) / (_1 * _2)).print(); printf("\n");
  std::cout << "Simplify test finished" << std::endl;
}
