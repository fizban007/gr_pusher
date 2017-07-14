#include "metrics.h"
#include "CudaLE.h"
#include "catch.hpp"

using namespace Aperture;
using CudaLE::D;
using CudaLE::simplify;

TEST_CASE("Printing metric name", "[metric]") {
  parse_metric("Cartesian");
  REQUIRE(metric::g_metric_cartesian.name() == "Cartesian");
}

TEST_CASE("Initialize metric", "[metric]") {
  // initialize_metric(MetricType::Schwarzschild, 4.0);
  parse_metric("Schwarzschild 4.0");
  REQUIRE(metric::g_metric_schwarzschild.m_rg == Approx(4.0));
  REQUIRE(metric::g_metric_schwarzschild.name() == "Schwarzschild4.000000");
}

TEST_CASE("log spherical", "[metric]") {
  parse_metric("Log_Spherical");
  double x1 = 2.0, x2 = 2.0;
  double r = std::exp(x1);
  double t = std::sin(x2);
  REQUIRE(metric::g_metric_log_spherical.g11(x1) - r * r == Approx(0.0));
  REQUIRE(metric::g_metric_log_spherical.g22(x1) - r * r == Approx(0.0));
  REQUIRE(metric::g_metric_log_spherical.g33(x1, x2) - r * r * t * t == Approx(0.0));

  // TODO: position to and from cartesian
}

TEST_CASE("Kerr Schild", "[metric]") {
  parse_metric("Kerr_Schild 1.0 0.9");
  REQUIRE(metric::g_metric_kerr_schild.m_rg == Approx(1.0));
  REQUIRE(metric::g_metric_kerr_schild.m_a == Approx(0.9));
}

TEST_CASE("is_GR", "[metric]") {
  REQUIRE(metric::g_metric_kerr_schild.is_GR() == true);
  REQUIRE(metric::g_metric_schwarzschild.is_GR() == true);
  REQUIRE(metric::g_metric_spherical.is_GR() == false);
  REQUIRE(metric::g_metric_log_spherical.is_GR() == false);
  REQUIRE(metric::g_metric_cartesian.is_GR() == false);

  metric::metric_base_interface* p = &metric::g_metric_kerr_schild;
  REQUIRE(p -> is_GR() == true);
}

TEST_CASE("Kerr Schild derivative", "[metric]") {
  parse_metric("Kerr_Schild 1.0 0.9");

  auto g = metric::g_metric_kerr_schild.g33;
  auto f = simplify(D<1>(D<2>(g)));
  auto f2 = simplify(D<2>(D<3>(g)));

  // f.print(); CudaLE::helper::print("\n");
  // f2.print(); CudaLE::helper::print("\n");

  REQUIRE(f(2.0, 1.3, 2.0) == Approx(2.06201));
  REQUIRE(f2(2.0, 1.3, 2.0) == Approx(0.0));
  // std::cout << f(2.0, 1.3, 2.0) << std::endl;
  // std::cout << f2(2.0, 1.3, 2.0) << std::endl;
}
