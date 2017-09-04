#include "solve.h"
#include <Eigen/Dense>
#include <iostream>

void solve(const std::array<std::array<double, 6>, 6>& J,
           const std::array<double, 6>& F, std::array<double, 6>& result) {
  Eigen::MatrixXd J_mat(6, 6);
  Eigen::VectorXd F_vec(6);
  for (int m = 0; m < 6; m++) {
    for (int n = 0; n < 6; n++) {
      J_mat(n, m) = J[n][m];
    }
    // F(m) = f[m].x();
    F_vec(m) = F[m];
  }
  // std::cout << J_mat << std::endl << std::endl;
  Eigen::VectorXd new_f = J_mat.colPivHouseholderQr().solve(F_vec);

  for (int i = 0; i < 6; i++) {
    result[i] = new_f(i);
  }
}

double norm(const std::array<double, 6>& f) {
  Eigen::VectorXd v(6);
  for (int i = 0; i < 6; i++) {
    v[i] = f[i];
  }
  return v.norm();
}
