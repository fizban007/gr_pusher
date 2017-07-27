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

TEST_CASE("Parsing a grid", "[grid]") {
  Grid g;

  g.parse({"DIM1 512 1.0 5.0 2",
          "DIM2 512 0.0 3.14 2",
          "DIM3 1   0.0 0.0  0"});

  REQUIRE(g.mesh().dims[0] == 516);
  REQUIRE(g.mesh().dims[1] == 516);
  REQUIRE(g.mesh().dims[2] == 1);
  REQUIRE(g.mesh().guard[0] == 2);
  REQUIRE(g.mesh().guard[1] == 2);
  REQUIRE(g.mesh().lower[0] == 1.0);
  REQUIRE(g.mesh().lower[1] == 0.0);

  REQUIRE(g.mesh().delta[0] == Approx(5.0 / 512));
  REQUIRE(g.mesh().delta[1] == Approx(3.14 / 512));
}

// TEST_CASE("Spherical face area", "[grid]") {
//   Grid g({"DIM1 64 1.0 5.00 2",
//           "DIM2 64 0.0 3.14 2",
//           "DIM3 1  0.0 6.28 0"});
//   // g.setup_areas(metrics::metric_spherical());
//   g.setup_areas(metric::g_metric_spherical, g, 30);

//   REQUIRE(g.face_area_array()[0].size() == 68 * 68);
//   REQUIRE(g.face_area(0, 30, 40, 0) - g.mesh().pos(0, 30, 1) * g.mesh().pos(0, 30, 1) *
//           (std::cos(g.mesh().pos(1, 39, 1)) - std::cos(g.mesh().pos(1, 40, 1))) * 6.28 == Approx(0.0).margin(1.0e-6));
//   REQUIRE(g.face_area(1, 30, 40, 0) -
//               (std::pow(g.mesh().pos(0, 30, 1), 3) - std::pow(g.mesh().pos(0, 29, 1), 3)) / 3.0 *
//           std::sin(g.mesh().pos(1, 40, 1)) * 6.28 == Approx(0.0).margin(1.0e-6));
//   REQUIRE(g.face_area(2, 30, 40, 0) - g.mesh().pos(0, 30, 0) *
//               0.5 * (g.mesh().pos(0, 30, 1) * g.mesh().pos(0, 30, 1) -
//                      g.mesh().pos(0, 29, 1) * g.mesh().pos(0, 29, 1)) *
//           std::sin(g.mesh().pos(1, 40, 0)) * g.mesh().delta[1] == Approx(0.0).margin(2.0e-6));

//   for (int i = 1; i < 67; i++) {
//     REQUIRE(g.face_area(1, i, 1, 0) == Approx(0.0).margin(1.0e-10));
//   }
// }

TEST_CASE("Scales in spherical coordinates", "[grid]") {
  Grid g({"DIM1 64 1.0 5.00 2",
          "DIM2 64 0.0 3.14 2",
          "DIM3 1  0.0 6.28 0"});
  g.setup_scales(metric::g_metric_spherical, g);

  REQUIRE(g.scales_array()[0].size() == 4);
  REQUIRE(g.scales_array()[0][0].size() == 68 * 68);
  for (int i = 2; i < 66; i++) {
    for (int j = 2; j < 66; j++) {
      double theta = g.mesh().pos(1, j, 1);
      double r = g.mesh().pos(0, i, 0);
      REQUIRE(g.scales(0, Stagger_t("001"))(i, j) == Approx(1.0));
      REQUIRE(g.scales(1, Stagger_t("010"))(i, j) == Approx(r));
      REQUIRE(g.scales(2, Stagger_t("010"))(i, j) == Approx(r * std::sin(theta)));
    }
  }
}

TEST_CASE("Scales in log spherical coordinates", "[grid]") {
  Grid g({"DIM1 64 0.0 2.00 2",
          "DIM2 64 0.0 3.14 2",
          "DIM3 1  0.0 6.28 0"});
  g.setup_scales(metric::g_metric_log_spherical, g);

  REQUIRE(g.scales_array()[0].size() == 4);
  REQUIRE(g.scales_array()[0][0].size() == 68 * 68);
  for (int i = 2; i < 66; i++) {
    for (int j = 2; j < 66; j++) {
      double r = g.mesh().pos(0, i, 0);
      double theta = g.mesh().pos(1, j, 1);
      REQUIRE(g.scales(0, Stagger_t("010"))(i, j) == Approx(std::exp(r)));
      REQUIRE(g.scales(1, Stagger_t("010"))(i, j) == Approx(std::exp(r)));
      REQUIRE(g.scales(2, Stagger_t("010"))(i, j) == Approx(std::exp(r) * std::sin(theta)));
    }
  }
}

TEST_CASE("Alpha and beta", "[grid]") {
  double a = 0.7, rg = 2.4;
  Grid g({"DIM1 64 1.0 5.00 2",
          "DIM2 64 0.0 3.14 2",
          "DIM3 1  0.0 0.01 0"});
  auto m = metric::metric_kerr_schild(rg, a);
  g.setup_metric(m, g);
  // g.setup_connection(m);
  // Eigen::Matrix4d g_uv;

  auto _rho2 = _r * _r + a * a * cos(_theta) * cos(_theta);
  auto _z = rg * _r / _rho2;
  auto _alpha = 1.0 / sqrt(1.0 + _z);

  VectorField<Scalar> beta_expected(g);
  // VectorField<Scalar> beta2_expected(g);
  ScalarField<Scalar> alpha_expected(g);
  beta_expected.initialize(0, _z / (1.0 + _z));
  // beta2_expected.initialize(1, _z / (1.0 + _z));
  alpha_expected.initialize(_alpha);
  // alpha_expected.initialize(1, _alpha);
  // alpha_expected.initialize(2, _alpha);
  for (int j = 1; j < 67; j++) {
    for (int i = 1; i < 67; i++) {
      Approx target = Approx(0.0).margin(1.0e-5);
      REQUIRE(g.beta_mask_array()[0] == 1);
      std::cout << "(" << i << ", " << j << ")" << std::endl;
      REQUIRE(g.beta_array()[0](i, j) == Approx(beta_expected(0, i, j)).margin(1.0e-5));
      // REQUIRE(g.beta2(1, i, j) - beta2_expected(1, i, j) == target);
      REQUIRE(g.alpha(i, j) - alpha_expected(i, j) == target);
      // REQUIRE(g.alpha(1, i, j) - alpha_expected(1, i, j) == target);
      // REQUIRE(g.alpha(2, i, j) - alpha_expected(2, i, j) == target);
    }
  }

  REQUIRE(g.beta_mask(0) == true);
  // REQUIRE(g.beta2_mask(2) == true);
  // REQUIRE(g.beta1_mask(0) == false);
  // REQUIRE(g.beta2_mask(1) == false);
  // REQUIRE(g.beta1_mask(2) == false);
  // REQUIRE(g.beta2_mask(0) == false);
}

TEST_CASE("Spherical alpha", "[grid]") {
  Grid g({"DIM1 64 1.0 5.00 2",
          "DIM2 64 0.0 3.14 2",
          "DIM3 1  0.0 0.01 0"});
  auto m = metric::metric_spherical();
  g.setup_metric(m, g);

  for (int j = 1; j < 67; j++) {
    for (int i = 1; i < 67; i++) {
      REQUIRE(g.alpha(i, j) == Approx(1.0));
      REQUIRE(g.alpha(i + j * 68, Vec3<double>(0.2 * g.mesh().delta[0],0.3 * g.mesh().delta[1],0.0)) == Approx(1.0));
    }
  }
}

TEST_CASE("Spherical connection", "[grid]") {

}

TEST_CASE("Making dual", "[grid]") {
  Grid g({"DIM1 32 0.0 2.00 2",
          "DIM2 32 0.0 3.14 2",
          "DIM3 32 0.0 6.28 2"});

  g.setup_metric(metric::g_metric_log_spherical, g);
  Grid g2 = g.make_dual();
  g2.setup_metric(metric::g_metric_log_spherical, g2);

  REQUIRE(g2.mesh().lower[0] - (g.mesh().lower[0] + g.mesh().delta[0] * 0.5) == Approx(0.0));
  REQUIRE(g2.mesh().lower[1] - (g.mesh().lower[1] + g.mesh().delta[1] * 0.5) == Approx(0.0));

  for (int i = 2; i < 34; i++) {
    for (int j = 2; j < 34; j++) {
      double r = g2.mesh().pos(0, i, 0);
      double r_s = g2.mesh().pos(0, i, 1);
      double theta = g2.mesh().pos(1, j, 0);
      double theta_s = g2.mesh().pos(1, j, 1);
      REQUIRE(g2.metric(0, 0, i, j) == Approx(std::exp(r) * std::exp(r)));
      REQUIRE(g2.metric(1, 1, i, j) == Approx(std::exp(r) * std::exp(r)));
      REQUIRE(g2.metric(2, 2, i, j) == Approx(std::exp(r) * std::sin(theta) * std::exp(r) * std::sin(theta)));
    }
  }
}

TEST_CASE("Saving and reloading", "[grid]") {
  Grid g({"DIM1 64 1.0 5.00 2",
          "DIM2 64 0.0 3.14 2",
          "DIM3 64 0.0 6.28 2"});
  g.setup_metric(metric::metric_spherical(), g);

  Grid g2({"DIM1 64 1.0 5.00 2",
          "DIM2 64 0.0 3.14 2",
          "DIM3 64 0.0 6.28 2"});
  g2.setup_metric(metric::metric_spherical(), g2);
  size_t N = g2.mesh().size();

  for (int n = 0; n < N; n++) {
    REQUIRE(g.det_array()[n] == g2.det_array()[n]);
    REQUIRE(g.alpha_array()[n] == g2.alpha_array()[n]);
  }
  for (int i = 0; i < 3; i++) {
    REQUIRE(g.mesh().lower[i] == g2.mesh().lower[i]);
    REQUIRE(g.mesh().sizes[i] == g2.mesh().sizes[i]);
    REQUIRE(g.mesh().delta[i] == g2.mesh().delta[i]);
    REQUIRE(g.mesh().guard[i] == g2.mesh().guard[i]);
    REQUIRE(g.mesh().dims[i] == g2.mesh().dims[i]);
    REQUIRE(g.mesh().dimension == g2.mesh().dimension);
    if (g.beta_mask_array()[i] == 1) {
      for (int n = 0; n < N; n++) {
        REQUIRE(g.beta_array()[i][n] == g2.beta_array()[i][n]);
        // REQUIRE(g.beta2_array()[i][n] == g2.beta2_array()[i][n]);
      }
    }
    for (int j = 0; j < 3; j++) {
      REQUIRE(g.metric_mask_array()[i][j] == g2.metric_mask_array()[i][j]);
      if (g.metric_mask_array()[i][j] == 1) {
        for (int n = 0; n < N; n++) {
          REQUIRE(g.metric_array()[i][j][n] == g2.metric_array()[i][j][n]);
        }
      }
    }
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        REQUIRE(g.conn_mask_array()[i][j][k] == g2.conn_mask_array()[i][j][k]);
        for (int n = 0; n < N; n++) {
          if (g.conn_mask_array()[i][j][k] == 1)
            REQUIRE(g.conn_array()[i][j][k][n] == g2.conn_array()[i][j][k][n]);
        }
      }
    }
  }
}
