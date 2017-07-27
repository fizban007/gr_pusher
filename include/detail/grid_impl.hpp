#ifndef _GRID_IMPL_HPP_
#define _GRID_IMPL_HPP_

#include <cctype>
#include <iostream>
#include <iomanip>
#include <utility>
#include <Eigen/Dense>
// #include <Eigen/Dense>
#include "grid.h"
#include "constant_defs.h"
#include "algorithms/interpolation.h"
#include "algorithms/interp_template.h"
#include "CudaLE.h"


// TODO: double check what kind of staggered definition for cell face area,
// determinant, metric elements makes the most sense, and can be made
// backward-compatible

#define SUM_PARTIAL_G(METRIC,A,I,MU,G_INV,POS)            \
  CudaLE::D<I + 1>( METRIC.g ## A ## 0 )(POS[0], POS[1], POS[2]) * G_INV[0][MU] \
  + CudaLE::D<I + 1>( METRIC.g ## A ## 1 )(POS[0], POS[1], POS[2]) * G_INV[1][MU] \
  + CudaLE::D<I + 1>( METRIC.g ## A ## 2 )(POS[0], POS[1], POS[2]) * G_INV[2][MU] \
  + CudaLE::D<I + 1>( METRIC.g ## A ## 3 )(POS[0], POS[1], POS[2]) * G_INV[3][MU]

#define SUM_DG(METRIC,I,MU,NU,G_INV,POS) \
  CudaLE::D<I + 1>( METRIC.g ## 0 ## 0 )(POS[0], POS[1], POS[2]) * G_INV[0][MU] * G_INV[0][NU] \
  + CudaLE::D<I + 1>( METRIC.g ## 0 ## 1 )(POS[0], POS[1], POS[2]) * G_INV[0][MU] * G_INV[1][NU] \
  + CudaLE::D<I + 1>( METRIC.g ## 0 ## 2 )(POS[0], POS[1], POS[2]) * G_INV[0][MU] * G_INV[2][NU] \
  + CudaLE::D<I + 1>( METRIC.g ## 0 ## 3 )(POS[0], POS[1], POS[2]) * G_INV[0][MU] * G_INV[3][NU] \
  + CudaLE::D<I + 1>( METRIC.g ## 1 ## 0 )(POS[0], POS[1], POS[2]) * G_INV[1][MU] * G_INV[0][NU] \
  + CudaLE::D<I + 1>( METRIC.g ## 1 ## 1 )(POS[0], POS[1], POS[2]) * G_INV[1][MU] * G_INV[1][NU] \
  + CudaLE::D<I + 1>( METRIC.g ## 1 ## 2 )(POS[0], POS[1], POS[2]) * G_INV[1][MU] * G_INV[2][NU] \
  + CudaLE::D<I + 1>( METRIC.g ## 1 ## 3 )(POS[0], POS[1], POS[2]) * G_INV[1][MU] * G_INV[3][NU] \
  + CudaLE::D<I + 1>( METRIC.g ## 2 ## 0 )(POS[0], POS[1], POS[2]) * G_INV[2][MU] * G_INV[0][NU] \
  + CudaLE::D<I + 1>( METRIC.g ## 2 ## 1 )(POS[0], POS[1], POS[2]) * G_INV[2][MU] * G_INV[1][NU] \
  + CudaLE::D<I + 1>( METRIC.g ## 2 ## 2 )(POS[0], POS[1], POS[2]) * G_INV[2][MU] * G_INV[2][NU] \
  + CudaLE::D<I + 1>( METRIC.g ## 2 ## 3 )(POS[0], POS[1], POS[2]) * G_INV[2][MU] * G_INV[3][NU] \
  + CudaLE::D<I + 1>( METRIC.g ## 3 ## 0 )(POS[0], POS[1], POS[2]) * G_INV[3][MU] * G_INV[0][NU] \
  + CudaLE::D<I + 1>( METRIC.g ## 3 ## 1 )(POS[0], POS[1], POS[2]) * G_INV[3][MU] * G_INV[1][NU] \
  + CudaLE::D<I + 1>( METRIC.g ## 3 ## 2 )(POS[0], POS[1], POS[2]) * G_INV[3][MU] * G_INV[2][NU] \
  + CudaLE::D<I + 1>( METRIC.g ## 3 ## 3 )(POS[0], POS[1], POS[2]) * G_INV[3][MU] * G_INV[3][NU]

namespace Aperture {

// const int N_subdivide = 1;
// const int N_subdivide = 30;

// Grid::Grid() {}
// Grid::~Grid() {}

// Grid::Grid(int N1, int N2, int N3) : m_mesh(N1, N2, N3), m_grid_hash(0),
//                                      m_metric_mask{}, m_connection_mask{} {}

// Grid::Grid(const std::vector<std::string>& config) :
//     m_metric_mask{},
//     m_connection_mask{} { parse(config); }

// Grid::Grid(const Grid& g) {
//   m_mesh = g.m_mesh;
//   m_grid_hash = g.m_grid_hash;
//   m_cache_filename = g.m_cache_filename;
//   for (int i = 0; i < 3; i++) {
//     m_cell_face_area[i] = g.m_cell_face_area[i];
//     m_cell_edge_length[i] = g.m_cell_edge_length[i];
//     // m_alpha[i] = g.m_alpha[i];
//     // m_beta1[i] = g.m_beta1[i];
//     // m_beta2[i] = g.m_beta2[i];
//     // m_det[i] = g.m_det[i];
//     for (unsigned int j = 0; j < 3; j++) {
//       m_metric[i][j] = g.m_metric[i][j];
//       m_metric_mask[i][j] = g.m_metric_mask[i][j];
//       // for (unsigned int k = 0; k < 3; k++) {
//       //   m_connection[i][j][k] = g.m_connection[i][j][k];
//       //   m_connection_mask[i][j][k] = g.m_connection_mask[i][j][k];
//       // }
//     }
//   }
// }

// Grid::Grid(Grid&& g)
//     : m_alpha(std::move(g.m_alpha)),
//       m_cell_face_area(std::move(g.m_cell_face_area)),
//       m_beta1(std::move(g.m_beta1)),
//       m_beta2(std::move(g.m_beta2)),
//       // m_cell_edge_length(std::move(g.m_cell_edge_length)),
//       m_det(std::move(g.m_det)),
//       m_metric(std::move(g.m_metric)),
//       m_metric_mask(std::move(g.m_metric_mask)),
//       m_connection(std::move(g.m_connection)),
//       m_connection_mask(std::move(g.m_connection_mask)) {
//   m_mesh = g.m_mesh;
//   m_grid_hash = g.m_grid_hash;
//   m_cache_filename = g.m_cache_filename;
//   // std::cout << "In move constructor!" << std::endl;
// }

// Grid& Grid::operator=(const Grid& g) {
//   m_mesh = g.m_mesh;
//   m_grid_hash = g.m_grid_hash;
//   m_cache_filename = g.m_cache_filename;
//   for (int i = 0; i < 3; i++) {
//     m_cell_face_area[i] = g.m_cell_face_area[i];
//     // m_cell_edge_length[i] = g.m_cell_edge_length[i];
//     m_alpha[i] = g.m_alpha[i];
//     m_beta1[i] = g.m_beta1[i];
//     m_beta2[i] = g.m_beta2[i];
//     m_det[i] = g.m_det[i];
//     for (unsigned int j = 0; j < 3; j++) {
//       m_metric[i][j] = g.m_metric[i][j];
//       m_metric_mask[i][j] = g.m_metric_mask[i][j];
//       for (unsigned int k = 0; k < 3; k++) {
//         m_connection[i][j][k] = g.m_connection[i][j][k];
//         m_connection_mask[i][j][k] = g.m_connection_mask[i][j][k];
//       }
//     }
//   }
//   return *this;
// }

// Grid& Grid::operator=(Grid&& g) {
//   m_mesh = g.m_mesh;
//   m_grid_hash = g.m_grid_hash;
//   m_cache_filename = g.m_cache_filename;
//   m_cell_face_area = std::move(g.m_cell_face_area);
//   // m_cell_edge_length = std::move(g.m_cell_edge_length);
//   m_alpha = std::move(g.m_alpha);
//   m_beta1 = std::move(g.m_beta1);
//   m_beta2 = std::move(g.m_beta2);
//   m_metric = std::move(g.m_metric);
//   m_metric_mask = std::move(g.m_metric_mask);
//   m_connection = std::move(g.m_connection);
//   m_connection_mask = std::move(g.m_connection_mask);
//   m_det = std::move(g.m_det);
//   // std::cout << "In move assignment!" << std::endl;
//   return *this;
// }

// template <typename Metric>
// void
// Grid::setup_lengths_f::operator()(const Metric& g) {
//   // Loop over the whole grid
//   double pos[3];
//   for (int k = 0; k < m_mesh.dims[2]; k++) {
//     pos[2] = m_mesh.pos(2, k, 1);
//     for (int j = 0; j < m_mesh.dims[1]; j++) {
//       pos[1] = m_mesh.pos(1, j, 1);
//       for (int i = 0; i < m_mesh.dims[0]; i++) {
//         pos[0] = m_mesh.pos(0, i, 1);
//         int idx = m_mesh.getIdx(i, j, k);

//         // First integrate side lengths
//         for (int n = 0; n < 3; n++) {
//           if (n < m_mesh.dim()) {
//             double length = 0.0;
//             pos[n] =
//                 pos[n] - m_mesh.delta[n] + m_mesh.delta[n] * 0.5 / N_subdivide;
//             for (int substep = 0; substep < N_subdivide; substep++) {
//               // length += sqrt(g(n + 1, n + 1, pos)) * m_mesh.delta[n] /
//               // N_subdivide;
//               length += g(n + 1, n + 1, pos) * m_mesh.delta[n] / N_subdivide;
//               pos[n] += m_mesh.delta[n] / N_subdivide;
//             }
//             pos[n] -= m_mesh.delta[n] * 0.5 / N_subdivide;
//             // m_cell_edge_length[n][idx] = length;
//           } else {
//             // m_cell_edge_length[n][idx] = sqrt(g(n + 1, n + 1, pos)) *
//             // m_mesh.sizes[n];
//             // m_cell_edge_length[n][idx] = g(n + 1, n + 1, pos) * m_mesh.sizes[n];
//           }
//         }
//       }
//     }
//   }
// }

// template <typename Metric>
// void
// Grid::setup_areas_f::operator()(const Metric& g, Grid& grid, int subdivide) const {
//   for (int i = 0; i < 3; i++) {
//     grid.m_cell_face_area[i].resize(grid.m_mesh.extent());
//   }
//   // Loop over the whole grid
//   double pos[3];
//   for (int k = 0; k < grid.m_mesh.dims[2]; k++) {
//     pos[2] = grid.m_mesh.pos(2, k, 1);
//     for (int j = 0; j < grid.m_mesh.dims[1]; j++) {
//       pos[1] = grid.m_mesh.pos(1, j, 1);
//       for (int i = 0; i < grid.m_mesh.dims[0]; i++) {
//         pos[0] = grid.m_mesh.pos(0, i, 1);
//         int idx = grid.m_mesh.get_idx(i, j, k);
//         // Then integrate face area
//         // FIXME: Is there any potential bugs in this part?
//         for (int n = 0; n < 3; n++) {
//           double area = 0.0;
//           int n_trans[2]{(n + 1) % 3, (n + 2) % 3};
//           if (n_trans[0] < grid.m_mesh.dim())
//             pos[n_trans[0]] = pos[n_trans[0]] - grid.m_mesh.delta[n_trans[0]] +
//                               grid.m_mesh.delta[n_trans[0]] * 0.5 / subdivide;
//           for (int step0 = 0; step0 < subdivide; step0++) {
//             if (n_trans[1] < grid.m_mesh.dim())
//               pos[n_trans[1]] = pos[n_trans[1]] - grid.m_mesh.delta[n_trans[1]] +
//                                 grid.m_mesh.delta[n_trans[1]] * 0.5 / subdivide;
//             for (int step1 = 0; step1 < subdivide; step1++) {
//               // +1 is to translate from 0 based to 1 based indexing
//               // area += g.det(pos) / std::sqrt(g(n + 1, n + 1, pos));
//               area += std::sqrt(g.det(pos));
//               if (n_trans[1] < grid.m_mesh.dim()) {
//                 pos[n_trans[1]] += grid.m_mesh.delta[n_trans[1]] / subdivide;
//               }
//             }
//             if (n_trans[0] < grid.m_mesh.dim())
//               pos[n_trans[0]] += grid.m_mesh.delta[n_trans[0]] / subdivide;
//             if (n_trans[1] < grid.m_mesh.dim())
//               pos[n_trans[1]] -= 0.5 * grid.m_mesh.delta[n_trans[1]] / subdivide;
//           }
//           if (n_trans[0] < grid.m_mesh.dim())
//             pos[n_trans[0]] -= 0.5 * grid.m_mesh.delta[n_trans[0]] / subdivide;
//           area *= grid.m_mesh.delta[n_trans[0]] / subdivide;
//           area *= grid.m_mesh.delta[n_trans[1]] / subdivide;
//           grid.m_cell_face_area[n][idx] = area;
//         }
//       }
//     }
//   }
// }

// template <typename Metric>
// void
// Grid::setup_alpha_beta(const Metric& g) {
//   // Loop over the whole grid
//   double pos[3];
//   for (int k = 0; k < m_mesh.dims[2]; k++) {
//     pos[2] = m_mesh.pos(2, k, 0);
//     for (int j = 0; j < m_mesh.dims[1]; j++) {
//       pos[1] = m_mesh.pos(1, j, 0);
//       for (int i = 0; i < m_mesh.dims[0]; i++) {
//         pos[0] = m_mesh.pos(0, i, 0);
//         int idx = m_mesh.getIdx(i, j, k);
//         // alpha are on face centers
//         for (int n = 0; n < 3; n++) {
//           int n_trans[2] = {(n + 1) % 3, (n + 2) % 3};
//           if (n < m_mesh.dim()) {
//             pos[n] += 0.5 * m_mesh.delta[n];
//             m_alpha[n][idx] = g.alpha(pos);
//             // double divisor[2] = {};
//             // FIXME: 1.0e-6 is a magic number!!!
//             // divisor[0] = std::max(std::sqrt(g(n_trans[1] + 1, n_trans[1] + 1,
//             // pos)), 1.0e-8);
//             // divisor[1] = std::max(std::sqrt(g(n_trans[0] + 1, n_trans[0] + 1,
//             // pos)), 1.0e-8);
//             // m_beta1[n][idx] = g.beta(n_trans[0] + 1, pos) * std::sqrt(g(n +
//             // 1, n + 1, pos)) / divisor[0];
//             // m_beta2[n][idx] = g.beta(n_trans[1] + 1, pos) * std::sqrt(g(n +
//             // 1, n + 1, pos)) / divisor[1];
//             m_beta1[n][idx] = g.beta(n_trans[0] + 1, pos);
//             m_beta2[n][idx] = g.beta(n_trans[1] + 1, pos);
//             // m_norm[n][idx] = std::sqrt(g(n + 1, n + 1, pos));
//             // m_norm[n][idx] = std::max(std::sqrt(g(n + 1, n + 1, pos)), 1.0e-6);
//             m_det[n][idx] = g.det(pos);
//             pos[n] -= 0.5 * m_mesh.delta[n];
//           } else {
//             m_alpha[n][idx] = g.alpha(pos);
//             // double divisor[2] = {};
//             // divisor[0] = std::max(std::sqrt(g(n_trans[1] + 1, n_trans[1] + 1,
//             // pos)), 1.0e-6);
//             // divisor[1] = std::max(std::sqrt(g(n_trans[0] + 1, n_trans[0] + 1,
//             // pos)), 1.0e-6);
//             // m_beta1[n][idx] = g.beta(n_trans[0] + 1, pos) * std::sqrt(g(n +
//             // 1, n + 1, pos)) / divisor[0];
//             // m_beta2[n][idx] = g.beta(n_trans[1] + 1, pos) * std::sqrt(g(n +
//             // 1, n + 1, pos)) / divisor[1];
//             m_beta1[n][idx] = g.beta(n_trans[0] + 1, pos);
//             m_beta2[n][idx] = g.beta(n_trans[1] + 1, pos);
//             // m_norm[n][idx] = std::sqrt(g(n + 1, n + 1, pos));
//             // m_norm[n][idx] = std::max(std::sqrt(g(n + 1, n + 1, pos)), 1.0e-6);
//             m_det[n][idx] = g.det(pos);
//           }
//         }
//       }
//     }
//   }
// }

// template <typename Metric>
// void
// Grid::setup_gamma(const Metric& g) {
//   // Loop over the whole grid
//   double pos[3];
//   for (int k = 0; k < m_mesh.dims[2]; k++) {
//     pos[2] = m_mesh.pos(2, k, 0);
//     for (int j = 0; j < m_mesh.dims[1]; j++) {
//       pos[1] = m_mesh.pos(1, j, 0);
//       for (int i = 0; i < m_mesh.dims[0]; i++) {
//         pos[0] = m_mesh.pos(0, i, 0);
//         int idx = m_mesh.getIdx(i, j, k);
//         // field components are on face centers
//         for (int n = 0; n < 3; n++) {
//           int n_trans[2] = {(n + 1) % 3, (n + 2) % 3};
//           if (n < m_mesh.dim()) {
//             pos[n] += 0.5 * m_mesh.delta[n];
//             for (int m = 0; m < 3; m++) {
//               if (m_metric[n][m].size() > 0) {
//                 m_metric[n][m][idx] = g(n + 1, m + 1, pos);
//               }
//             }
//             pos[n] -= 0.5 * m_mesh.delta[n];
//           } else {
//             for (int m = 0; m < 3; m++) {
//               if (m_metric[n][m].size() > 0)
//                 m_metric[n][m][idx] = g(n + 1, m + 1, pos);
//             }
//           }
//         }
//       }
//     }
//   }
// }

template <typename Metric>
void
Grid::setup_metric_f::operator()(const Metric& g, Grid& grid) const {
  // grid.m_type = (MetricType)Metric::type;
  grid.m_metric_ptr = &g;
  grid.m_type = grid.m_metric_ptr -> type();
  // Need to adjust the grid hash according to the metric
  std::size_t metric_hash = grid.hash_line(g.name());
  grid.m_grid_hash <<= 1;
  grid.m_grid_hash ^= metric_hash;
  grid.gen_file_name();

  // Allocate det, alpha, and beta arrays
  grid.allocate_arrays();
  // Mask the corresponding beta array
  // if (g.b1 != CudaLE::ZeroOp())
  //   grid.m_beta1_mask[1] = grid.m_beta2_mask[2] = true;
  // if (g.b2 != CudaLE::ZeroOp())
  //   grid.m_beta1_mask[2] = grid.m_beta2_mask[0] = true;
  // if (g.b3 != CudaLE::ZeroOp())
  //   grid.m_beta1_mask[0] = grid.m_beta2_mask[1] = true;

  // Resize the metric arrays, always do the diagonal
  grid.m_metric[0][0].resize(grid.m_mesh.extent());
  grid.m_metric[1][1].resize(grid.m_mesh.extent());
  grid.m_metric[2][2].resize(grid.m_mesh.extent());
  grid.m_metric_mask[0][0] = grid.m_metric_mask[1][1] = grid.m_metric_mask[2][2] = 1;
  // optionally do the non-diagonal depending on spatial metric
  if (g.g12 != CudaLE::ZeroOp()) {
    grid.m_metric[0][1].resize(grid.m_mesh.extent());
    grid.m_metric[1][0].resize(grid.m_mesh.extent());
    grid.m_metric_mask[0][1] = grid.m_metric_mask[1][0] = 1;
  }
  if (g.g13 != CudaLE::ZeroOp()) {
    grid.m_metric[0][2].resize(grid.m_mesh.extent());
    grid.m_metric[2][0].resize(grid.m_mesh.extent());
    grid.m_metric_mask[0][2] = grid.m_metric_mask[2][0] = 1;
  }
  if (g.g23 != CudaLE::ZeroOp()) {
    grid.m_metric[1][2].resize(grid.m_mesh.extent());
    grid.m_metric[2][1].resize(grid.m_mesh.extent());
    grid.m_metric_mask[1][2] = grid.m_metric_mask[2][1] = 1;
  }

  if (g.b1 != CudaLE::ZeroOp()) grid.m_beta_mask[0] = 1;
  if (g.b2 != CudaLE::ZeroOp()) grid.m_beta_mask[1] = 1;
  if (g.b3 != CudaLE::ZeroOp()) grid.m_beta_mask[2] = 1;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      grid.m_inv_metric[i][j].resize(grid.m_mesh.extent());
    }
  }

  // Resize the connection arrays
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        grid.m_connection[i][j][k].resize(grid.m_mesh.dims[0], grid.m_mesh.dims[1], grid.m_mesh.dims[2]);
      }
    }
  }

  std::ifstream f_test;
  f_test.open(grid.m_cache_filename.c_str());
  if (f_test.good()) {
    f_test.close();
    grid.load_from_disk();
  } else {
    f_test.close();
    // setup_areas(g);
    // setup_gamma(g);
    // setup_connection(g);

    // Loop over the whole grid
    double pos[3];
    for (int k = 0; k < grid.m_mesh.dims[2]; k++) {
      pos[2] = grid.m_mesh.pos(2, k, 0);
      for (int j = 0; j < grid.m_mesh.dims[1]; j++) {
        pos[1] = grid.m_mesh.pos(1, j, 0);
        for (int i = 0; i < grid.m_mesh.dims[0]; i++) {
          pos[0] = grid.m_mesh.pos(0, i, 0);
          int idx = grid.m_mesh.get_idx(i, j, k);
          // Eigen::Matrix3d g_mat = Eigen::Matrix3d::Zero();
          for (int n = 0; n < 3; n++) {
            // int n_trans[2] = {(n + 1) % 3, (n + 2) % 3};
            if (n < grid.m_mesh.dim()) {
              // We define metric coefficients at face centers with respect to the
              // first index
              for (int m = 0; m < 3; m++) {
                if (grid.m_metric[n][m].size() > 0) {
                  grid.m_metric[n][m][idx] = g(n + 1, m + 1, pos);
                }
              }
              // Determinant is also defined at cell faces
              // pos[n] += 0.5 * grid.m_mesh.delta[n];
              grid.m_det[idx] = sqrt(g.det(pos));
              grid.m_alpha[idx] = g.alpha(pos);
              if (grid.m_beta_mask[n] == 1)
                grid.m_beta[n][idx] = g.beta(n, pos);
              // grid.m_beta2[n][idx] = g.beta(n_trans[1] + 1, pos);
              // pos[n] -= 0.5 * grid.m_mesh.delta[n];
            } else {
              for (int m = 0; m < 3; m++) {
                // metric coefficients are defined at the grid centers
                if (grid.m_metric[n][m].size() > 0) {
                  grid.m_metric[n][m][idx] = g(n + 1, m + 1, pos);
                  // g_mat(n, m) = grid.m_metric[n][m][idx];
                }
              }
              grid.m_det[idx] = sqrt(g.det(pos));
              grid.m_alpha[idx] = g.alpha(pos);
              if (grid.m_beta_mask[n] == 1)
                grid.m_beta[n][idx] = g.beta(n, pos);
              // grid.m_beta2[n][idx] = g.beta(n_trans[1] + 1, pos);
            }
          }
          // Eigen::Matrix3d g_inv = g_mat.inverse();
          // for (int n = 0; n < 3; n++) {
          //   for (int m = 0; m < 3; m++) {
          //     grid.m_inv_metric[n][m][idx] = g_inv(n, m);
          //   }
          // }
          // }
        }
      }
    }
    grid.setup_inv_metric();
    grid.setup_connection(g);

    grid.save_to_disk();
  }
  // std::cout << "Min resolved length is " << min_resolved_length() << std::endl;
}

template <typename Metric>
void
Grid::setup_scales_f::operator()(const Metric& g, Grid& grid) const {
  auto mesh = grid.m_mesh;
  int dim = mesh.dim();
  int stagger_num = (1 << dim);
  for (int i = 0; i < 3; i++) {
    grid.m_scales[i].resize(stagger_num);
    for (int j = 0; j < stagger_num; j++) {
      grid.m_scales[i][j].resize(mesh.extent());
    }
  }

  for (int n = 0; n < stagger_num; n++) {
    // Using bit operations to be more intuitive
    int istag = check_bit(n, 0); // Extract lowest bit
    int jstag = check_bit(n, 1); // Extract second lowest bit
    int kstag = check_bit(n, 2); // Extract highest bit
    // loop over cell indices(i,j,k)
    for (int k = 0; k < mesh.dims[2]; k++) {
      for (int j = 0; j < mesh.dims[1]; j++) {
        for (int i = 0; i < mesh.dims[0]; i++) {
          // calculate the coordinate values
          double q1 = mesh.pos(0, i, istag) + EPS;
          double q2 = mesh.pos(1, j, jstag) + EPS;
          double q3 = mesh.pos(2, k, kstag) + EPS;
          // calculate the scaling functions h1, h2, h3
          grid.m_scales[0][n](i, j, k) = std::sqrt(g.g11(q1, q2, q3));
          grid.m_scales[1][n](i, j, k) = std::sqrt(g.g22(q1, q2, q3));
          grid.m_scales[2][n](i, j, k) = std::sqrt(g.g33(q1, q2, q3));
        }
      }
    }

  }
}

template <typename Metric>
void
Grid::setup_connection(const Metric& g) {
  // Loop over the whole grid
  double pos[3];
  // Eigen::Matrix4d g_uv, g_inv;
  std::array<std::array<double, 4>, 4> g_uv, g_inv;
  for (int k = 0; k < m_mesh.dims[2]; k++) {
    pos[2] = m_mesh.pos(2, k, 0);
    for (int j = 0; j < m_mesh.dims[1]; j++) {
      pos[1] = m_mesh.pos(1, j, 0);
      for (int i = 0; i < m_mesh.dims[0]; i++) {
        pos[0] = m_mesh.pos(0, i, 0);
        // connection components are at cell centers
        // we first construct the matrix g_{\mu\nu}
        for (int n = 0; n < 4; n++) {
          for (int m = 0; m < 4; m++) {
            g_uv[n][m] = g(n, m, pos);
            // g_0m += g_uv(n, m) * g.beta(m, pos);
          }
        }
        // g_inv = g_uv.inverse();
        // invert_matrix(g_uv, g_inv);
        g_inv = invert_mat4(g_uv);
        // if (i == 20 && j == 20 && k == 0) {
        //   std::cout << g_inv << std::endl;
        // }

        int idx = m_mesh.get_idx(i, j, k);
        for (int mu = 0; mu < 4; mu++) {
          for (int nu = 0; nu < 4; nu++) {
            m_connection[0][mu][nu][idx] = SUM_DG(g,0,mu,nu,g_inv,pos);
            if (m_connection[0][mu][nu][idx] != 0.0) m_connection_mask[0][mu][nu] = 1;
            m_connection[1][mu][nu][idx] = SUM_DG(g,1,mu,nu,g_inv,pos);
            if (m_connection[1][mu][nu][idx] != 0.0) m_connection_mask[1][mu][nu] = 1;
            m_connection[2][mu][nu][idx] = SUM_DG(g,2,mu,nu,g_inv,pos);
            if (m_connection[2][mu][nu][idx] != 0.0) m_connection_mask[2][mu][nu] = 1;
          }
          // This is the only loop we can iterate over
          // First index is i, second is alpha, in the expression
          // g_{\alpha\beta, i}g^{\beta\mu}
          // m_connection[0][0][mu] = PARTIAL_G(g,0,0,1,pos) * g_inv(0, mu);
        //   m_connection[0][0][mu][idx] = SUM_DG(g,0,0,mu,g_inv,pos);
        //   if (m_connection[0][0][mu][idx] != 0.0) m_connection_mask[0][0][mu] = 1;
        //   m_connection[0][1][mu][idx] = SUM_DG(g,1,0,mu,g_inv,pos);
        //   if (m_connection[0][1][mu][idx] != 0.0) m_connection_mask[0][1][mu] = 1;
        //   m_connection[0][2][mu][idx] = SUM_DG(g,2,0,mu,g_inv,pos);
        //   if (m_connection[0][2][mu][idx] != 0.0) m_connection_mask[0][2][mu] = 1;
        //   m_connection[0][3][mu][idx] = SUM_DG(g,3,0,mu,g_inv,pos);
        //   if (m_connection[0][3][mu][idx] != 0.0) m_connection_mask[0][3][mu] = 1;
        //   m_connection[1][0][mu][idx] = SUM_DG(g,0,1,mu,g_inv,pos);
        //   if (m_connection[1][0][mu][idx] != 0.0) m_connection_mask[1][0][mu] = 1;
        //   m_connection[1][1][mu][idx] = SUM_DG(g,1,1,mu,g_inv,pos);
        //   if (m_connection[1][1][mu][idx] != 0.0) m_connection_mask[1][1][mu] = 1;
        //   m_connection[1][2][mu][idx] = SUM_DG(g,2,1,mu,g_inv,pos);
        //   if (m_connection[1][2][mu][idx] != 0.0) m_connection_mask[1][2][mu] = 1;
        //   m_connection[1][3][mu][idx] = SUM_DG(g,3,1,mu,g_inv,pos);
        //   if (m_connection[1][3][mu][idx] != 0.0) m_connection_mask[1][3][mu] = 1;
        //   // By default the simulation is 2D, therefore third derivative is
        //   // always zero because we assume it's symmetry direction
        //   if (m_mesh.dim() > 2) {
        //     m_connection[2][0][mu][idx] = SUM_DG(g,0,2,mu,g_inv,pos);
        //     if (m_connection[2][0][mu][idx] != 0.0) m_connection_mask[2][0][mu] = 1;
        //     m_connection[2][1][mu][idx] = SUM_DG(g,1,2,mu,g_inv,pos);
        //     if (m_connection[2][1][mu][idx] != 0.0) m_connection_mask[2][1][mu] = 1;
        //     m_connection[2][2][mu][idx] = SUM_DG(g,2,2,mu,g_inv,pos);
        //     if (m_connection[2][2][mu][idx] != 0.0) m_connection_mask[2][2][mu] = 1;
        //     m_connection[2][3][mu][idx] = SUM_DG(g,3,2,mu,g_inv,pos);
        //     if (m_connection[2][3][mu][idx] != 0.0) m_connection_mask[2][3][mu] = 1;
        //   }
        }
      }
    }
  }
}

template <typename Metric>
Grid
Grid::make_dual_f::operator()(const Metric& g, bool inverse) const {
  Grid grid = grid.make_dual(inverse);
  grid.setup_metric(g, grid);
  return grid;
}


// template <int InterpolationOrder, typename POS_T>
// Matrix3
// Grid::metric_matrix(int cell, const Vec3<POS_T>& rel_pos) const {
//   interpolator<InterpolationOrder> interp;
//   Matrix3 result = Matrix3::ZeroOp();

//   Vec3<int> c = m_mesh.getCell3D(cell);
//   Vec3<int> lower = c - interp.radius();
//   Vec3<int> upper = c + interp.support() - interp.radius();
//   if (dim() < 3) {
//     lower[2] = upper[2] = c[2];
//   }
//   for (int n = 0; n < 3; n++) {
//     for (int m = 0; m < 3; m++) {
//       if (m_metric_mask[n][m] != 1) continue;
//       for (int k = lower[2]; k <= upper[2]; k++) {
//         for (int j = lower[1]; j <= upper[1]; j++) {
//           for (int i = lower[0]; i <= upper[0]; i++) {
//             if (dim() < 3) {
//               // if (m_metric[n][m].size() > 0) {
//                 result(n, m) += m_metric[n][m](i, j, k) * interp.interp_cell(rel_pos[0], c[0], i, (n == 0 ? 1 : 0))
//                                 * interp.interp_cell(rel_pos[1], c[1], j, (n == 1 ? 1 : 0));
//               // }
//             } else {
//               result(n, m) += m_metric[n][m](i, j, k) * interp.interp_cell(rel_pos[0], c[0], i, (n == 0 ? 1 : 0))
//                               * interp.interp_cell(rel_pos[1], c[1], j, (n == 1 ? 1 : 0))
//                               * interp.interp_cell(rel_pos[2], c[2], k, (n == 2 ? 1 : 0));
//             }
//           }
//         }
//       }
//     }
//   }
//   return result;
// }

// template <int InterpolationOrder, typename POS_T>
// Matrix3
// Grid::metric_matrix(const Vec3<POS_T>& pos) const {
//   Vec3<POS_T> rel_pos;
//   int cell = m_mesh.findCell(pos, rel_pos);
//   return metric_matrix<InterpolationOrder>(cell, rel_pos);
// }

// template <typename POS_T>
// Scalar
// Grid::alpha(int cell, const Vec3<POS_T>& rel_pos) const {
//   Interpolator interp(1);
//   // Original metric array has stagger same as vector component n
//   Scalar result = 0.0;

//   Vec3<int> c = m_mesh.get_cell_3d(cell);
//   Vec3<int> lower = c - interp.radius();
//   Vec3<int> upper = c + interp.support() - interp.radius();
//   if (dim() < 3) {
//     lower[2] = upper[2] = c[2];
//   }
//   for (int k = lower[2]; k <= upper[2]; k++) {
//     for (int j = lower[1]; j <= upper[1]; j++) {
//       for (int i = lower[0]; i <= upper[0]; i++) {
//         result += m_alpha[0](i, j, k) * interp.interp_cell(rel_pos[0], c[0], i, 1)
//                   * interp.interp_cell(rel_pos[1], c[1], j, 0)
//                   * (dim() < 3 ? 1.0 : interp.interp_cell(rel_pos[2], c[2], k, 0));
//       }
//     }
//   }
//   return result;
// }

template <typename POS_T>
Scalar
Grid::alpha(const Vec3<POS_T>& pos) const {
  Vec3<POS_T> rel_pos;
  int cell = m_mesh.find_cell(pos, rel_pos);
  return alpha(cell, rel_pos);
}

// template <int InterpolationOrder, typename POS_T>
// Scalar
// Grid::det(int cell, const Vec3<POS_T>& rel_pos) const {
//   interpolator<InterpolationOrder> interp;
//   // Original metric array has stagger same as vector component n
//   Scalar result = 0.0;

//   Vec3<int> c = m_mesh.getCell3D(cell);
//   Vec3<int> lower = c - interp.radius();
//   Vec3<int> upper = c + interp.support() - interp.radius();
//   if (dim() < 3) {
//     lower[2] = upper[2] = c[2];
//   }
//   for (int k = lower[2]; k <= upper[2]; k++) {
//     for (int j = lower[1]; j <= upper[1]; j++) {
//       for (int i = lower[0]; i <= upper[0]; i++) {
//         result += m_det[0](i, j, k) * interp.interp_cell(rel_pos[0], c[0], i, 1)
//                   * interp.interp_cell(rel_pos[1], c[1], j, 0)
//                   * (dim() < 3 ? 1.0 : interp.interp_cell(rel_pos[2], c[2], k, 0));
//       }
//     }
//   }
//   return result;
// }

// template <int InterpolationOrder, typename POS_T>
// Scalar
// Grid::det(const Vec3<POS_T>& pos) const {
//   Vec3<POS_T> rel_pos;
//   int cell = m_mesh.findCell(pos, rel_pos);
//   return det<InterpolationOrder>(cell, rel_pos);
// }

// template <int InterpolationOrder, typename POS_T>
// Vector3
// Grid::beta(int cell, const Vec3<POS_T>& rel_pos) const {
//   interpolator<InterpolationOrder> interp;
//   // Original metric array has stagger same as vector component n
//   Vector3 result(0.0, 0.0, 0.0);

//   Vec3<int> c = m_mesh.getCell3D(cell);
//   Vec3<int> lower = c - interp.radius();
//   Vec3<int> upper = c + interp.support() - interp.radius();
//   if (dim() < 3) {
//     lower[2] = upper[2] = c[2];
//   }
//   for (int k = lower[2]; k <= upper[2]; k++) {
//     for (int j = lower[1]; j <= upper[1]; j++) {
//       for (int i = lower[0]; i <= upper[0]; i++) {
//         if (dim() < 3) {
//           result[0] += m_beta2[1](i, j, k) * interp.interp_cell(rel_pos[0], c[0], i, 0)
//                        * interp.interp_cell(rel_pos[1], c[1], j, 1);
//           result[1] += m_beta1[0](i, j, k) * interp.interp_cell(rel_pos[0], c[0], i, 1)
//                        * interp.interp_cell(rel_pos[1], c[1], j, 0);
//           result[2] += m_beta2[0](i, j, k) * interp.interp_cell(rel_pos[0], c[0], i, 1)
//                        * interp.interp_cell(rel_pos[1], c[1], j, 0);
//         } else {
//           result[0] += m_beta1[2](i, j, k) * interp.interp_cell(rel_pos[0], c[0], i, 0)
//                        * interp.interp_cell(rel_pos[1], c[1], j, 0)
//                        * interp.interp_cell(rel_pos[2], c[2], k, 1);
//           result[1] += m_beta1[0](i, j, k) * interp.interp_cell(rel_pos[0], c[0], i, 1)
//                        * interp.interp_cell(rel_pos[1], c[1], j, 0)
//                        * interp.interp_cell(rel_pos[2], c[2], k, 0);
//           result[2] += m_beta1[1](i, j, k) * interp.interp_cell(rel_pos[0], c[0], i, 0)
//                        * interp.interp_cell(rel_pos[1], c[1], j, 1)
//                        * interp.interp_cell(rel_pos[2], c[2], k, 0);
//         }
//       }
//     }
//   }
//   return result;
// }

// template <int InterpolationOrder, typename POS_T>
// Vector3
// Grid::beta(const Vec3<POS_T>& pos) const {
//   Vec3<POS_T> rel_pos;
//   int cell = m_mesh.findCell(pos, rel_pos);
//   return beta<InterpolationOrder>(cell, rel_pos);
// }

// template <int InterpolationOrder, typename POS_T>
// void
// Grid::connection(int cell, const Vec3<POS_T>& rel_pos, double conn[3][4][4]) const {
//   interpolator<InterpolationOrder> interp;
//   // All connections are evaluated at the center of the grid

//   Vec3<int> c = m_mesh.getCell3D(cell);
//   Vec3<int> lower = c - interp.radius();
//   Vec3<int> upper = c + interp.support() - interp.radius();
//   if (dim() < 3) {
//     lower[2] = upper[2] = c[2];
//   }
//   for (int a = 0; a < 3; a++) {
//     for (int m = 0; m < 4; m++) {
//       for (int n = 0; n < 4; n++) {
//         conn[a][m][n] = 0.0;
//       }
//     }
//   }
//   for (int a = 0; a < 3; a++) {
//     for (int m = 0; m < 4; m++) {
//       for (int n = 0; n < 4; n++) {
//         if (m_connection_mask[a][m][n] != 1) continue;
//         for (int k = lower[2]; k <= upper[2]; k++) {
//           for (int j = lower[1]; j <= upper[1]; j++) {
//             for (int i = lower[0]; i <= upper[0]; i++) {
//               conn[a][m][n] += m_connection[a][m][n](i, j, k) * interp.interp_cell(rel_pos[0], c[0], i, 0)
//                                * interp.interp_cell(rel_pos[1], c[1], j, 0)
//                                * (dim() < 3 ? 1.0 : interp.interp_cell(rel_pos[2], c[2], k, 0));
//             }
//           }
//         }
//       }
//     }
//   }
// }

// template <int InterpolationOrder, typename POS_T>
// void
// Grid::connection(const Vec3<POS_T>& pos, double conn[3][4][4]) const {
//   Vec3<POS_T> rel_pos;
//   int cell = m_mesh.findCell(pos, rel_pos);
//   connection<InterpolationOrder>(cell, rel_pos, conn);
// }

// Scalar
// Grid::min_resolved_length(int cell) const {
//   // Scalar lambda_p = std::sqrt(m_metric[0][0][cell]) * m_mesh.delta[0];
//   // if (dim() > 1) {
//   //   Scalar tmp =  std::sqrt(m_metric[1][1][cell]) * m_mesh.delta[1];
//   //   lambda_p = lambda_p > tmp ? lambda_p : tmp;
//   // }
//   // if (dim() > 2) {
//   //   Scalar tmp =  std::sqrt(m_metric[2][2][cell]) * m_mesh.delta[2];
//   //   lambda_p = lambda_p > tmp ? lambda_p : tmp;
//   // }
//   double lambda_p = std::min(m_mesh.delta[0], m_mesh.delta[1]);
//   return lambda_p;
// }

// Scalar
// Grid::cell_volume(int cell) const {
//   return m_det[2][cell] * m_mesh.delta[0] * m_mesh.delta[1] * (dim() > 2 ? m_mesh.delta[2] : 1.0);
// }

template <typename Double>
Double
Grid::det(int cell, const Double& x1, const Double& x2, const Double& x3) const {
  Double result = 0.0;

  Vec3<int> c = m_mesh.get_cell_3d(cell);
  Vec3<int> lower = c - Vec3<int>(1, 1, 1);
  Vec3<int> upper = c + Vec3<int>(1, 1, 1);
  if (dim() < 3) {
    lower[2] = upper[2] = c[2];
  }
  for (int k = lower[2]; k <= upper[2]; k++) {
    for (int j = lower[1]; j <= upper[1]; j++) {
      for (int i = lower[0]; i <= upper[0]; i++) {
        result += m_det(i, j, k) * interp_cell(x1, i - c[0], m_mesh.delta[0], 0)
                  * interp_cell(x2, j - c[1], m_mesh.delta[1], 0)
                  * (dim() < 3 ? 1.0 : interp_cell(x3, k - c[2], m_mesh.delta[2], 0));
      }
    }
  }
  return result;
}

template <typename Double>
Double
Grid::alpha(int cell, const Double& x1, const Double& x2, const Double& x3) const {
  Double result = 0.0;

  Vec3<int> c = m_mesh.get_cell_3d(cell);
  Vec3<int> lower = c - Vec3<int>(1, 1, 1);
  Vec3<int> upper = c + Vec3<int>(1, 1, 1);
  if (dim() < 3) {
    lower[2] = upper[2] = c[2];
  }
  for (int k = lower[2]; k <= upper[2]; k++) {
    for (int j = lower[1]; j <= upper[1]; j++) {
      for (int i = lower[0]; i <= upper[0]; i++) {
        result += m_alpha(i, j, k) * interp_cell(x1, i - c[0], m_mesh.delta[0], 0)
                  * interp_cell(x2, j - c[1], m_mesh.delta[1], 0)
                  * (dim() < 3 ? 1.0 : interp_cell(x3, k - c[2], m_mesh.delta[2], 0));
      }
    }
  }
  return result;
}

template <typename Double>
Double
Grid::beta(int n, int cell, const Double& x1, const Double& x2, const Double& x3) const {
  Double result = 0.0;

  Vec3<int> c = m_mesh.get_cell_3d(cell);
  Vec3<int> lower = c - Vec3<int>(1, 1, 1);
  Vec3<int> upper = c + Vec3<int>(1, 1, 1);
  if (dim() < 3) {
    lower[2] = upper[2] = c[2];
  }
  for (int k = lower[2]; k <= upper[2]; k++) {
    for (int j = lower[1]; j <= upper[1]; j++) {
      for (int i = lower[0]; i <= upper[0]; i++) {
        if (dim() < 3) {
          // if (n == 0)
          result += m_beta[n](i, j, k) * interp_cell(x1, i - c[0], m_mesh.delta[0], 0)
                      * interp_cell(x2, j - c[1], m_mesh.delta[1], 0);
          // else if (n == 1)
          //   result += m_beta1[0](i, j, k) * interp_cell(x1, i - c[0], m_mesh.delta[0], 1)
          //             * interp_cell(x2, j - c[1], m_mesh.delta[1], 0);
          // else if (n == 2)
          //   result += m_beta2[0](i, j, k) * interp_cell(x1, i - c[0], m_mesh.delta[0], 1)
          //             * interp_cell(x2, j - c[1], m_mesh.delta[1], 0);
        } else {
          // if (n == 0)
          result += m_beta[n](i, j, k) * interp_cell(x1, i - c[0], m_mesh.delta[0], 0)
                    * interp_cell(x2, j - c[1], m_mesh.delta[1], 0)
                    * interp_cell(x3, k - c[2], m_mesh.delta[2], 0);
                       // * interp.interp_cell(rel_pos[2], c[2], k, 1);
          // else if (n == 1)
          //   result += m_beta1[0](i, j, k) * interp_cell(x1, i - c[0], m_mesh.delta[0], 1)
          //             * interp_cell(x2, j - c[1], m_mesh.delta[1], 0)
          //             * interp_cell(x3, k - c[2], m_mesh.delta[2], 0);
          //              // * interp.interp_cell(rel_pos[1], c[1], j, 0)
          //              // * interp.interp_cell(rel_pos[2], c[2], k, 0);
          // else if (n == 2)
          //   result += m_beta1[1](i, j, k) * interp_cell(x1, i - c[0], m_mesh.delta[0], 0)
          //             * interp_cell(x2, j - c[1], m_mesh.delta[1], 1)
          //             * interp_cell(x3, k - c[2], m_mesh.delta[2], 0);
          //              // * interp.interp_cell(rel_pos[1], c[1], j, 1)
          //              // * interp.interp_cell(rel_pos[2], c[2], k, 0);
        }
      }
    }
  }
  return result;
}

template <typename Double>
Double
Grid::metric(int n, int m, int cell, const Double& x1, const Double& x2, const Double& x3) const {
  Double result = 0.0;

  Vec3<int> c = m_mesh.get_cell_3d(cell);
  Vec3<int> lower = c - Vec3<int>(1, 1, 1);
  Vec3<int> upper = c + Vec3<int>(1, 1, 1);
  if (dim() < 3) {
    lower[2] = upper[2] = c[2];
  }
  for (int k = lower[2]; k <= upper[2]; k++) {
    for (int j = lower[1]; j <= upper[1]; j++) {
      for (int i = lower[0]; i <= upper[0]; i++) {
        if (m_metric_mask[n][m] == 1)
          result += m_metric[n][m](i, j, k) * interp_cell(x1, i - c[0], m_mesh.delta[0], 0)
                    * interp_cell(x2, j - c[1], m_mesh.delta[1], 0)
                    * (dim() < 3 ? 1.0 : interp_cell(x3, k - c[2], m_mesh.delta[2], 0));
      }
    }
  }
  return result;
}

template <typename Double>
Double
Grid::inv_metric(int n, int m, int cell, const Double& x1, const Double& x2, const Double& x3) const {
  Double result = 0.0;

  Vec3<int> c = m_mesh.get_cell_3d(cell);
  Vec3<int> lower = c - Vec3<int>(1, 1, 1);
  Vec3<int> upper = c + Vec3<int>(1, 1, 1);
  if (dim() < 3) {
    lower[2] = upper[2] = c[2];
  }
  for (int k = lower[2]; k <= upper[2]; k++) {
    for (int j = lower[1]; j <= upper[1]; j++) {
      for (int i = lower[0]; i <= upper[0]; i++) {
        result += m_inv_metric[n][m](i, j, k) * interp_cell(x1, i - c[0], m_mesh.delta[0], 0)
                  * interp_cell(x2, j - c[1], m_mesh.delta[1], 0)
                  * (dim() < 3 ? 1.0 : interp_cell(x3, k - c[2], m_mesh.delta[2], 0));
      }
    }
  }
  return result;
}

template <typename Double>
Double
Grid::connection(int n, int u, int v, int cell, const Double& x1, const Double& x2, const Double& x3) const {
  Double result = 0.0;

  Vec3<int> c = m_mesh.get_cell_3d(cell);
  Vec3<int> lower = c - Vec3<int>(1, 1, 1);
  Vec3<int> upper = c + Vec3<int>(1, 1, 1);
  if (dim() < 3) {
    lower[2] = upper[2] = c[2];
  }
  for (int k = lower[2]; k <= upper[2]; k++) {
    for (int j = lower[1]; j <= upper[1]; j++) {
      for (int i = lower[0]; i <= upper[0]; i++) {
        if (m_connection_mask[n][u][v] == 1)
          result += m_connection[n][u][v](i, j, k) * interp_cell(x1, i - c[0], m_mesh.delta[0], 0)
                    * interp_cell(x2, j - c[1], m_mesh.delta[1], 0)
                    * (dim() < 3 ? 1.0 : interp_cell(x3, k - c[2], m_mesh.delta[2], 0));
      }
    }
  }
  return result;
}

template <typename Double>
Double
Grid::det(int cell, const Vec3<Double>& x) const {
  return det(cell, x[0], x[1], x[2]);
}

template <typename Double>
Double
Grid::alpha(int cell, const Vec3<Double>& x) const {
  return alpha(cell, x[0], x[1], x[2]);
}

template <typename Double>
Double
Grid::beta(int n, int cell, const Vec3<Double>& x) const {
  return beta(n, cell, x[0], x[1], x[2]);
}

template <typename Double>
Double
Grid::metric(int i, int j, int cell, const Vec3<Double>& x) const {
  return metric(i, j, cell, x[0], x[1], x[2]);
}

template <typename Double>
Double
Grid::inv_metric(int i, int j, int cell, const Vec3<Double>& x) const {
  return inv_metric(i, j, cell, x[0], x[1], x[2]);
}

template <typename Double>
Double
Grid::connection(int i, int u, int v, int cell, const Vec3<Double>& x) const {
  return connection(i, u, v, cell, x[0], x[1], x[2]);
}

}

#ifdef SUM_PARTIAL_G
#undef SUM_PARTIAL_G
#endif

#ifdef SUM_DG
#undef SUM_DG
#endif

#endif  // _GRID_IMPL_HPP_
