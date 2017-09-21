#include "catch.hpp"
#include "CudaLE.h"
#include <random>
#include <vector>
#include <iostream>
#include "vec3.h"

// using namespace Aperture;
using namespace CudaLE;
using namespace CudaLE::placeholders;

std::default_random_engine generator;
std::uniform_real_distribution<float> dist(0.0, 1.0);

template <typename T1, typename T2>
bool random_tests(const T1& f1, const T2& f2, int num_samples = 1000, int args = 1) {
  // std::vector<float> nums(num_samples);
  bool result = true;
  for (int i = 0; i < num_samples; i++) {
    float x = dist(generator);
    float y1 = f1(x), y2 = f2(x);
    if (!(y1 == Approx(y2))) {
      std::cout << "At x = " << x << ", f1 != f2!" << std::endl;
      std::cout << "y1 = " << y1 << ", y2 = " << y2 << std::endl;
      return false;
    }
  }
  return true;
}

template <typename T1, typename T2>
bool random_tests3(const T1& f1, const T2& f2, int num_samples = 1000) {
  // std::vector<float> nums(num_samples);
  bool result = true;
  for (int i = 0; i < num_samples; i++) {
    Aperture::Vec3<float> x;
    x[0] = dist(generator);
    x[1] = dist(generator);
    x[2] = dist(generator);
    float y1 = f1(x[0], x[1], x[2]), y2 = f2(x[0], x[1], x[2]);
    if (!(y1 == Approx(y2))) {
      std::cout << "At x = " << x << ", f1 != f2!" << std::endl;
      std::cout << "y1 = " << y1 << ", y2 = " << y2 << std::endl;
      return false;
    }
  }
  return true;
}



TEST_CASE("Testing OneOp and ZeroOp", "[CudaLE]") {
  auto f = ConstOp(3.0);
  REQUIRE(D<1>(f) == ZeroOp());

  auto g = placeholders::_2;
  REQUIRE(D<1>(g) == ZeroOp());
  REQUIRE(D<2>(g) == OneOp());

  REQUIRE(D<1>(4.0) == ZeroOp());

  REQUIRE(D<1>(ZeroOp()) == ZeroOp());

}

TEST_CASE("Testing Algebraic Binary Ops", "[CudaLE]") {
  // Add
  auto f1 = _1 + square(_1);
  REQUIRE(random_tests(D<1>(f1), OneOp{} + 2.0 * _1));
  REQUIRE(random_tests(D<2>(f1), ZeroOp{}));

  // Minus
  auto f2 = _1 - cos(_1);
  REQUIRE(random_tests(D<1>(f2), OneOp{} + sin(_1)));

  // Multiply
  auto f3 = _1 * exp(_1);
  REQUIRE(random_tests(D<1>(f3), exp(_1) + _1 * exp(_1)));

  // Divide
  auto f4 = exp(_1) / _1;
  REQUIRE(random_tests(D<1>(f4), (exp(_1) * _1 - exp(_1)) / square(_1)));
}

TEST_CASE("Testing Unary Ops", "[CudaLE]") {
  auto f1 = sin(_1);
  REQUIRE(random_tests(D<1>(f1), [](float x) {return cos(x);}));

  auto f2 = cos(_1);
  REQUIRE(random_tests(D<1>(f2), [](float x) {return -sin(x);}));

  auto f3 = exp(_1);
  REQUIRE(random_tests(D<1>(f3), [](float x) {return exp(x);}));

  auto f4 = sqrt(_1);
  REQUIRE(random_tests(D<1>(f4), [](float x) {return 0.5 / sqrt(x);}));

  auto f5 = square(_1);
  REQUIRE(random_tests(D<1>(f5), [](float x) {return 2.0 * x;}));

  auto f6 = log(_1);
  REQUIRE(random_tests(D<1>(f6), [](float x) {return 1.0 / x;}));

  auto f7 = pow<3>(_1);
  REQUIRE(random_tests(D<1>(f7), [](float x) {return 3.0 * x * x;}));

  auto f8 = ConstOp(3.0) * sin(_1);
  REQUIRE(random_tests(D<1>(f8), [](float x) {return 3.0 * cos(x);}));
}

TEST_CASE("Composite operators", "[CudaLE]") {
  auto f1 = sin(exp(_1));
  REQUIRE(random_tests(D<1>(f1), exp(_1) * cos(exp(_1))));

  auto f2 = log(square(exp(_1)));
  REQUIRE(random_tests(D<1>(f2), 2.0 * exp(_1) * exp(_1) / square(exp(_1))));

  typedef double Data;
  CudaLE::Var<1, Data> _x;
  CudaLE::Var<2, Data> _y;
  CudaLE::Var<3, Data> _z;
  double a_ = 0.65;
  DEFINE_FUNCTOR(R2, square(_x) + square(_y) + square(_z));
  DEFINE_FUNCTOR(r,
                 sqrt(0.5 * (R2 - a_ * a_ + sqrt(square(R2 - a_ * a_) +
                                                 4.0 * a_ * a_ * square(_z)))));

  REQUIRE(random_tests3(D<1>(r), (2.0 * _x + ((2.0 * _x * (-a_ * a_ + R2)) / sqrt(
            square(2.0 * a_ * _z) + square(-a_ * a_ + R2)))) /
                       (sqrt(8.0) * sqrt(-a_ * a_ + R2 + sqrt(square(2.0 * a_ * _z) + square(-a_ * a_ + R2))))));
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
