#ifndef _GRID_H_
#define _GRID_H_

#include <array>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
// #include <functional>
#include "constant_defs.h"
#include "multi_array.h"
#include "quadmesh.h"
#include "typedefs.h"
#include "vec3.h"
#include "metrics.h"

namespace Aperture {

class Grid {
 public:
  Grid();
  Grid(int N1, int N2 = 1, int N3 = 1);
  Grid(const std::array<std::string, 3>& config);
  Grid(const Grid& g);
  Grid(Grid&& g);
  ~Grid();

  void parse(const std::array<std::string, 3>& config);

  Grid& operator=(const Grid& g);
  Grid& operator=(Grid&& g);
  bool operator==(const Grid& g) const { return (m_mesh == g.m_mesh); }

  struct setup_metric_f {
    template <typename Metric>
    void operator()(const Metric& g, Grid& grid) const;
  } setup_metric;

  // struct setup_lenghts_f {
  //   template <typename Metric>
  //   void operator()(const Metric& g, Grid& grid, int subdivide = 1) const;
  // } setup_lengths;

  struct setup_areas_f {
    template <typename Metric>
    void operator()(const Metric& g, Grid& grid, int subdivide = 1) const;
  } setup_areas;

  struct setup_scales_f {
    template <typename Metric>
    void operator()(const Metric& g, Grid& grid) const;
  } setup_scales;

  void setup_inv_metric();

  // template <typename Metric>
  // void setup_alpha_beta(const Metric& g);

  // template <typename Metric>
  // void setup_gamma(const Metric& g);

  template <typename Metric>
  void setup_connection(const Metric& g);

  MultiArray<Scalar>& scales(int n, Stagger_t stagger);
  const MultiArray<Scalar>& scales(int n, Stagger_t stagger) const;

  Scalar metric(int n, int m, int cell) const { return m_metric[n][m][cell]; }
  Scalar metric(int n, int m, int c1, int c2, int c3 = 0) const {
    return m_metric[n][m](c1, c2, c3);
  }
  bool metric_mask(int n, int m) const { return m_metric_mask[n][m]; }
  bool beta1_mask(int n) const { return m_beta1_mask[n]; }
  bool beta2_mask(int n) const { return m_beta2_mask[n]; }
  bool conn_mask(int n, int u, int v) const { return m_connection_mask[n][u][v] == 1; }

  // Scalar metric(int n, int m, const Vec3<Scalar>& global_pos,
  //               int interp = 1) const;
  // Scalar metric(int n, int m, int cell, const Vec3<Pos_t>& rel_pos,
  //               int interp = 1) const;

  // Matrix3 metric_matrix(int cell, const Vec3<Pos_t>& rel_pos, int interp = 1)
  // const;
  // Matrix3 metric_matrix(const Vec3<Pos_t>& global_pos, int interp = 1) const;

  Scalar norm(int n, int cell) const { return std::sqrt(m_metric[n][n][cell]); }
  Scalar norm(int n, int c1, int c2 = 0, int c3 = 0) const {
    return std::sqrt(m_metric[n][n](c1, c2, c3));
  }

  Scalar inner_product(const Vec3<Pos_t>& pos, const Vec3<Scalar>& v1,
                       const Vec3<Scalar>& v2, int interp = 1) const;

  // template <typename Metric>
  // void setup_beta(const Metric& g);

  // Scalar face_area(int n, int cell) const { return m_cell_face_area[n][cell]; }
  // Scalar face_area(int n, int c1, int c2, int c3 = 0) const {
    // return m_cell_face_area[n](c1, c2, c3);
  // }
  // Scalar edge_length(int n, int cell) const {
  //   return m_cell_edge_length[n][cell];
  // }
  // Scalar edge_length(int n, int c1, int c2, int c3 = 0) const {
  //   return m_cell_edge_length[n](c1, c2, c3);
  // }

  Scalar min_resolved_length(int cell) const;
  Scalar cell_volume(int cell) const;

  Scalar alpha(int n, int cell) const { return m_alpha[n][cell]; }
  // template <typename POS_T>
  // Scalar alpha(int cell, const Vec3<POS_T>& rel_pos) const;
  template <typename POS_T>
  Scalar alpha(const Vec3<POS_T>& pos) const;

  Scalar beta1(int n, int cell) const { return m_beta1[n][cell]; }
  Scalar beta2(int n, int cell) const { return m_beta2[n][cell]; }
  // Vector3 beta(int cell, const Vec3<Pos_t>& rel_pos) const;
  // Vector3 beta(const Vec3<Scalar>& pos) const;

  // These methods are for interpolation and AD
  template <typename Double>
  Double det(int cell, const Double& x1, const Double& x2, const Double& x3) const;
  template <typename Double>
  Double alpha(int cell, const Double& x1, const Double& x2, const Double& x3) const;
  template <typename Double>
  Double beta(int n, int cell, const Double& x1, const Double& x2, const Double& x3) const;
  template <typename Double>
  Double metric(int i, int j, int cell, const Double& x1, const Double& x2, const Double& x3) const;
  template <typename Double>
  Double inv_metric(int i, int j, int cell, const Double& x1, const Double& x2, const Double& x3) const;
  template <typename Double>
  Double connection(int i, int u, int v, int cell, const Double& x1, const Double& x2, const Double& x3) const;
  template <typename Double>
  Double det(int cell, const Vec3<Double>& x) const;
  template <typename Double>
  Double alpha(int cell, const Vec3<Double>& x) const;
  template <typename Double>
  Double beta(int n, int cell, const Vec3<Double>& x) const;
  template <typename Double>
  Double metric(int i, int j, int cell, const Vec3<Double>& x) const;
  template <typename Double>
  Double inv_metric(int i, int j, int cell, const Vec3<Double>& x) const;
  template <typename Double>
  Double connection(int i, int u, int v, int cell, const Vec3<Double>& x) const;

  Scalar det(int n, int cell) const { return m_det[n][cell]; }
  // template <typename Pos_t>
  // Scalar det(int cell, const Vec3<Pos_t>& rel_pos) const;
  // template <typename Pos_t>
  // Scalar det(const Vec3<Pos_t>& pos) const;

  // First index is i, second is alpha, in the expression
  // g_{\alpha\beta, i}g^{\beta\mu}
  Scalar connection(int i, int a, int u, int pos) const {
    return m_connection[i][a][u][pos];
  }
  template <typename Pos_t>
  void connection(int cell, const Vec3<Pos_t>& rel_pos, double conn[3][4][4],
  int interp = 1) const;
  template <typename Pos_t>
  void connection(const Vec3<Pos_t>& pos, double conn[3][4][4], int interp =
  1) const;

  Scalar connection(int i, int u, int a, int c1, int c2, int c3 = 0) const {
    return m_connection[i][a][u](c1, c2, c3);
  }

  Scalar alpha(int n, int c1, int c2, int c3 = 0) const {
    return m_alpha[n](c1, c2, c3);
  }
  Scalar beta1(int n, int c1, int c2, int c3 = 0) const {
    return m_beta1[n](c1, c2, c3);
  }
  Scalar beta2(int n, int c1, int c2, int c3 = 0) const {
    return m_beta2[n](c1, c2, c3);
  }
  Scalar det(int n, int c1, int c2, int c3 = 0) const {
    return m_det[n](c1, c2, c3);
  }

  Grid make_dual(bool inverse = false) const;

  struct make_dual_f {
    template <typename Metric>
    Grid operator()(const Metric& g, bool from_dual = false) const;
  } make_dual_m;
  // Grid make_dual(const Metric& g, bool from_dual = false) const;

  std::array<std::string, 3> gen_config() const;

  Vec3<Pos_t> rel_pos_in_dual(const Vec3<Pos_t>& pos, int& cell) const;

  Quadmesh& mesh() { return m_mesh; }
  const Quadmesh& mesh() const { return m_mesh; }
  int size() const { return m_mesh.size(); }
  Extent extent() const { return m_mesh.extent(); }
  unsigned int dim() const { return m_mesh.dim(); }
  std::size_t grid_hash() const { return m_grid_hash; }
  MetricType type() const { return m_type; }
  bool isGR() const { return m_metric_ptr->is_GR(); }

  // const std::array<MultiArray<Scalar>, 3>& face_area_array() const {
  //   return m_cell_face_area;
  // }
  // const std::array<multi_array<Scalar>, 3>& edge_length_array() const {
  //   return m_cell_edge_length;
  // }
  // const std::array<multi_array<Scalar>, 3>& beta1_array() const {
  // return m_beta1;
  // }
  // const std::array<multi_array<Scalar>, 3>& beta2_array() const {
  //   return m_beta2;
  // }
  // const std::array<multi_array<Scalar>, 3>& alpha_array() const {
  //   return m_alpha;
  // }
  const std::array<MultiArray<Scalar>, 3>& det_array() const { return m_det; }
  const std::array<MultiArray<Scalar>, 3>& alpha_array() const { return m_alpha; }
  const std::array<MultiArray<Scalar>, 3>& beta1_array() const { return m_beta1; }
  const std::array<MultiArray<Scalar>, 3>& beta2_array() const { return m_beta2; }
  const std::array<std::array<char, 3>, 3>& metric_mask_array() const {
    return m_metric_mask;
  }
  const std::array<std::array<MultiArray<Scalar>, 3>, 3>& metric_array() const {
    return m_metric;
  }
  const std::array<std::vector<MultiArray<Scalar>>, 3>& scales_array() const {
    return m_scales;
  }
  const std::array<std::array<std::array<MultiArray<Scalar>, 4>, 4>, 3>&
  conn_array() const {
    return m_connection;
  }
  const std::array<std::array<std::array<char, 4>, 4>, 3>&
  conn_mask_array() const {
    return m_connection_mask;
  }

  void save_to_disk();
  void load_from_disk();
  void gen_file_name();

 private:
  void parse_line(const std::string& line);
  void parse_line_nums(std::stringstream& nums, int n);
  void allocate_arrays();
  std::array<std::array<double, 4>, 4> invert_mat4(const std::array<std::array<double, 4>, 4>& mat);

  std::size_t hash_line(const std::string& line, std::size_t seed = 0);
  std::size_t hash(const std::array<std::string, 3>& str, std::size_t seed = 0);

  Quadmesh m_mesh;
  std::size_t m_grid_hash;
  std::string m_cache_filename;
  MetricType m_type;
  const metric::metric_base_interface* m_metric_ptr;

  // alpha and beta are defined on face-centered positions
  // std::array<MultiArray<Scalar>, 3> m_cell_face_area;
  std::array<MultiArray<Scalar>, 3> m_alpha;
  // We only define two components of beta because we only need two transverse
  // components of beta at every cell face
  std::array<MultiArray<Scalar>, 3> m_beta1;
  std::array<MultiArray<Scalar>, 3> m_beta2;
  std::array<bool, 3> m_beta1_mask;
  std::array<bool, 3> m_beta2_mask;
  // std::array<multi_array<Scalar>, 3> m_cell_edge_length;
  std::array<MultiArray<Scalar>, 3> m_det;

  // TODO: work out the stagger of the metric tensor
  std::array<std::array<MultiArray<Scalar>, 3>, 3> m_metric;
  std::array<std::array<MultiArray<Scalar>, 3>, 3> m_inv_metric;
  std::array<std::array<char, 3>, 3> m_metric_mask;

  // When using the differential solver it is better to use scale functions
  std::array<std::vector<MultiArray<Scalar>>, 3> m_scales;

  // The outmost array has dimension 3 because we assume metric is
  // time-independent
  // These coefficients are defined by g_{\alpha\beta, i}g^{\beta\mu}
  // \alpha and \mu are free indices from 0 to 3, and i can only be 1 to 3.
  std::array<std::array<std::array<MultiArray<Scalar>, 4>, 4>, 3>
  m_connection;
  std::array<std::array<std::array<char, 4>, 4>, 3> m_connection_mask;
};
}

#include "detail/grid_impl.hpp"

#endif  // _GRID_H_
