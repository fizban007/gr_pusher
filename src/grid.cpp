#include "grid.h"
#include "CudaLE.h"
#include "algorithms/interpolation.h"
#include <Eigen/Dense>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <utility>

using namespace Aperture;

Grid::Grid() :
    Grid(1, 1, 1) {}

Grid::~Grid() {}

Grid::Grid(int N1, int N2, int N3)
    : m_mesh(N1, N2, N3),
      m_grid_hash(0),
      m_beta_mask{0, 0, 0},
      // m_beta2_mask{},
      m_metric_mask{} {
  // allocate_arrays();
}
// , m_connection_mask{} {}

Grid::Grid(const std::array<std::string, 3>& config)
    : m_beta_mask{}, m_metric_mask{} {
  parse(config);
}
// m_connection_mask{} { parse(config); }

Grid::Grid(const Grid& g) {
  m_mesh = g.m_mesh;
  m_grid_hash = g.m_grid_hash;
  m_cache_filename = g.m_cache_filename;
  m_alpha = g.m_alpha;
  m_det = g.m_det;
  for (int i = 0; i < 3; i++) {
    // m_cell_face_area[i] = g.m_cell_face_area[i];
    // m_cell_edge_length[i] = g.m_cell_edge_length[i];
    m_beta[i] = g.m_beta[i];
    m_beta_mask[i] = g.m_beta_mask[i];
    for (unsigned int j = 0; j < 3; j++) {
      m_metric[i][j] = g.m_metric[i][j];
      m_inv_metric[i][j] = g.m_inv_metric[i][j];
      m_metric_mask[i][j] = g.m_metric_mask[i][j];
      for (unsigned int k = 0; k < 3; k++) {
        m_connection[i][j][k] = g.m_connection[i][j][k];
        m_connection_mask[i][j][k] = g.m_connection_mask[i][j][k];
      }
    }
  }
}

Grid::Grid(Grid&& g)
    : // m_cell_face_area(std::move(g.m_cell_face_area)),
      m_alpha(std::move(g.m_alpha)),
      // m_beta1(std::move(g.m_beta1)),
      m_beta(std::move(g.m_beta)),
      m_beta_mask(std::move(g.m_beta_mask)),
      // m_beta2_mask(std::move(g.m_beta2_mask)),
      // m_cell_edge_length(std::move(g.m_cell_edge_length)),
      m_det(std::move(g.m_det)),
      m_metric(std::move(g.m_metric)),
      m_inv_metric(std::move(g.m_inv_metric)),
      m_metric_mask(std::move(g.m_metric_mask)),
      m_connection(std::move(g.m_connection)),
      m_connection_mask(std::move(g.m_connection_mask)) {
  m_mesh = g.m_mesh;
  m_grid_hash = g.m_grid_hash;
  m_cache_filename = g.m_cache_filename;
  // std::cout << "In move constructor!" << std::endl;
}

Grid&
Grid::operator=(const Grid& g) {
  m_mesh = g.m_mesh;
  m_grid_hash = g.m_grid_hash;
  m_cache_filename = g.m_cache_filename;
  m_alpha = g.m_alpha;
  m_det = g.m_det;
  for (int i = 0; i < 3; i++) {
    // m_cell_face_area[i] = g.m_cell_face_area[i];
    // m_cell_edge_length[i] = g.m_cell_edge_length[i];
    m_beta[i] = g.m_beta[i];
    m_beta_mask[i] = g.m_beta_mask[i];
    // m_beta2[i] = g.m_beta2[i];
    // m_beta2_mask[i] = g.m_beta2_mask[i];
    for (unsigned int j = 0; j < 3; j++) {
      m_metric[i][j] = g.m_metric[i][j];
      m_inv_metric[i][j] = g.m_inv_metric[i][j];
      m_metric_mask[i][j] = g.m_metric_mask[i][j];
      // for (unsigned int k = 0; k < 3; k++) {
      //   m_connection[i][j][k] = g.m_connection[i][j][k];
      //   m_connection_mask[i][j][k] = g.m_connection_mask[i][j][k];
      // }
    }
  }
  return *this;
}

Grid&
Grid::operator=(Grid&& g) {
  m_mesh = g.m_mesh;
  m_grid_hash = g.m_grid_hash;
  m_cache_filename = g.m_cache_filename;
  // m_cell_face_area = std::move(g.m_cell_face_area);
  // m_cell_edge_length = std::move(g.m_cell_edge_length);
  m_alpha = std::move(g.m_alpha);
  m_beta = std::move(g.m_beta);
  m_beta_mask = std::move(g.m_beta_mask);
  m_metric = std::move(g.m_metric);
  m_inv_metric = std::move(g.m_inv_metric);
  m_metric_mask = std::move(g.m_metric_mask);
  m_connection = std::move(g.m_connection);
  m_connection_mask = std::move(g.m_connection_mask);
  m_det = std::move(g.m_det);
  // std::cout << "In move assignment!" << std::endl;
  return *this;
}

void
Grid::setup_inv_metric() {
  for (int k = 0; k < m_mesh.dims[2]; k++) {
    for (int j = 0; j < m_mesh.dims[1]; j++) {
      for (int i = 0; i < m_mesh.dims[0]; i++) {
        int idx = m_mesh.get_idx(i, j, k);
        Eigen::Matrix3d g_mat = Eigen::Matrix3d::Zero();
        for (int n = 0; n < 3; n++) {
          for (int m = 0; m < 3; m++) {
            if (m_metric[n][m].size() > 0) {
              g_mat(n, m) = m_metric[n][m][idx];
            }
          }
        }
        Eigen::Matrix3d g_inv = g_mat.inverse();

        for (int n = 0; n < 3; n++) {
          for (int m = 0; m < 3; m++) {
            m_inv_metric[n][m][idx] = g_inv(n, m);
          }
        }
      }
    }
  }
}

Scalar
Grid::min_resolved_length(int cell) const {
  // Scalar lambda_p = std::sqrt(m_metric[0][0][cell]) * m_mesh.delta[0];
  // if (dim() > 1) {
  //   Scalar tmp =  std::sqrt(m_metric[1][1][cell]) * m_mesh.delta[1];
  //   lambda_p = lambda_p > tmp ? lambda_p : tmp;
  // }
  // if (dim() > 2) {
  //   Scalar tmp =  std::sqrt(m_metric[2][2][cell]) * m_mesh.delta[2];
  //   lambda_p = lambda_p > tmp ? lambda_p : tmp;
  // }
  double lambda_p = std::min(m_mesh.delta[0], m_mesh.delta[1]);
  return lambda_p;
}

// FIXME: This is almost always a good choice but could mess up when the third
// direction is not symmetric or not a killing vector. What do?
Scalar
Grid::cell_volume(int cell) const {
  return m_det[cell] * m_mesh.delta[0] * m_mesh.delta[1] *
         (dim() > 2 ? m_mesh.delta[2] : 1.0);
}

void
Grid::parse_line(const std::string& line) {
  std::stringstream parseline(line);
  std::string word;
  parseline >> word;
  std::transform(word.begin(), word.end(), word.begin(), ::toupper);

  if (word.compare("DIM1") == 0) {
    parse_line_nums(parseline, 0);
  } else if (word.compare("DIM2") == 0) {
    parse_line_nums(parseline, 1);
  } else if (word.compare("DIM3") == 0) {
    parse_line_nums(parseline, 2);
  }
}

void
Grid::parse_line_nums(std::stringstream& nums, int n) {
  std::string size_s, n_s, guard_s, lower_s;
  nums >> n_s >> lower_s >> size_s >> guard_s;
  if (size_s != "") m_mesh.sizes[n] = std::atof(size_s.c_str());
  if (guard_s != "") m_mesh.guard[n] = std::atoi(guard_s.c_str());
  if (lower_s != "") m_mesh.lower[n] = std::atof(lower_s.c_str());
  if (n_s != "") m_mesh.dims[n] = std::atoi(n_s.c_str()) + 2 * m_mesh.guard[n];
}

void
Grid::parse(const std::array<std::string, 3>& config) {
  for (auto const& line : config) {
    parse_line(line);
  }
  hash(config);
  gen_file_name();

  // Figure out grid spacing from mesh dimensions and sizes
  for (int i = 0; i < 3; i++) {
    if (m_mesh.dims[i] > 1) {
      m_mesh.delta[i] = m_mesh.sizes[i] / m_mesh.reduced_dim(i);
    } else {
      m_mesh.delta[i] = m_mesh.sizes[i];
    }
  }
  m_mesh.dimension = m_mesh.dim();

  // allocate_arrays();
}

Grid
Grid::make_dual(bool inverse) const {
  Grid g(*this);
  for (unsigned int i = 0; i < 3; i++) {
    if (i < g.dim()) {
      // Making dual grid from regular grid, shifting everything by half a grid
      if (!inverse) {
        g.m_mesh.lower[i] += 0.5 * g.m_mesh.delta[i];
      } else {
        g.m_mesh.lower[i] -= 0.5 * g.m_mesh.delta[i];
      }
    }
  }
  g.hash(g.gen_config());
  g.gen_file_name();
  return g;
}

std::size_t
Grid::hash_line(const std::string& line, std::size_t seed) {
  std::size_t result = seed;
  for (std::string::size_type i = 0; i < line.size(); i++) {
    result = result * 101 + line[i];
  }
  // std::cout << "Hashing " << line << std::endl;

  return result;
}

std::size_t
Grid::hash(const std::array<std::string, 3>& config, std::size_t seed) {
  m_grid_hash = 0;
  for (auto const& line : config) {
    m_grid_hash <<= 1;
    m_grid_hash ^= hash_line(line);
  }
  return m_grid_hash;
}

void
Grid::gen_file_name() {
  std::stringstream ss;
  ss << ".gridcache." << m_grid_hash;
  m_cache_filename = ss.str();
}

std::array<std::string, 3>
Grid::gen_config() const {
  std::array<std::string, 3> result;
  for (int i = 0; i < 3; i++) {
    std::stringstream ss;
    ss.precision(4);
    ss << "DIM" << i + 1 << " " << m_mesh.reduced_dim(i) << " "
       << m_mesh.lower[i] << " " << m_mesh.sizes[i] << " " << m_mesh.guard[i];
    result[i] = ss.str();
  }
  return result;
}

void
Grid::save_to_disk() {
  std::cout << "Saving cache to file " << m_cache_filename << std::endl;
  std::ofstream fs;
  fs.open(m_cache_filename.c_str(), std::ios::out | std::ios::binary);
  fs << m_mesh;
  fs.write((char*)m_alpha.data(), m_mesh.size() * sizeof(Scalar));
  fs.write((char*)m_det.data(), m_mesh.size() * sizeof(Scalar));
  for (int i = 0; i < 3; i++) {
    // fs.write((char*)m_cell_edge_length[i].data(),
    //          m_mesh.size() * sizeof(Scalar));
    // fs.write((char*)m_cell_face_area[i].data(), m_mesh.size() * sizeof(Scalar));
    fs.write(&m_beta_mask[i], sizeof(char));
    if (m_beta_mask[i] == 1)
      fs.write((char*)m_beta[i].data(), m_mesh.size() * sizeof(Scalar));
    // fs.write((char*)m_norm[i].data(), m_mesh.size() * sizeof(Scalar));
    for (int j = 0; j < 3; j++) {
      fs.write(&m_metric_mask[i][j], sizeof(char));
      if (m_metric_mask[i][j] == 1)
        fs.write((char*)m_metric[i][j].data(), m_mesh.size() * sizeof(Scalar));
    }
    for (int j = 0; j < 3; j++) {
      fs.write((char*)m_inv_metric[i][j].data(), m_mesh.size() * sizeof(Scalar));
    }
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        fs.write(&m_connection_mask[i][j][k], sizeof(char));
        if (m_connection_mask[i][j][k] == 1)
          fs.write((char*)m_connection[i][j][k].data(), m_mesh.size() * sizeof(Scalar));
        // fs.write((char*)m_connection_mask[i][j][k], m_mesh.size() * sizeof(char));
      }
    }
  }
  fs.close();
}

void
Grid::load_from_disk() {
  std::cout << "Loading cache from file " << m_cache_filename << std::endl;
  std::ifstream fs;
  fs.open(m_cache_filename.c_str(), std::ios::in | std::ios::binary);
  fs >> m_mesh;
  for (int i = 0; i < 3; i++) {
    if (m_mesh.dims[i] > 1) {
      m_mesh.delta[i] = m_mesh.sizes[i] / m_mesh.reduced_dim(i);
    } else {
      m_mesh.delta[i] = m_mesh.sizes[i];
    }
  }
  char c;
  fs.get(c);  // hack to get an additional position on fs
  // std::cout << (int)fs.tellg() << std::endl;
  fs.read((char*)m_alpha.data(), m_mesh.size() * sizeof(Scalar));
  fs.read((char*)m_det.data(), m_mesh.size() * sizeof(Scalar));
  for (int i = 0; i < 3; i++) {
    // fs.read((char*)m_cell_edge_length[i].data(),
    //         m_mesh.size() * sizeof(Scalar));
    // fs.read((char*)m_cell_face_area[i].data(), m_mesh.size() * sizeof(Scalar));
    // fs.read((char*)m_alpha[i].data(), m_mesh.size() * sizeof(Scalar));
    fs.read(&m_beta_mask[i], sizeof(char));
    if (m_beta_mask[i] == 1)
      fs.read((char*)m_beta[i].data(), m_mesh.size() * sizeof(Scalar));
    // fs.read((char*)m_beta2[i].data(), m_mesh.size() * sizeof(Scalar));
    // fs.read((char*)m_norm[i].data(), m_mesh.size() * sizeof(Scalar));
    // fs.read((char*)m_det[i].data(), m_mesh.size() * sizeof(Scalar));
    for (int j = 0; j < 3; j++) {
      fs.read(&m_metric_mask[i][j], sizeof(char));
      if (m_metric_mask[i][j] == 1)
        fs.read((char*)m_metric[i][j].data(), m_mesh.size() * sizeof(Scalar));
    }
    for (int j = 0; j < 3; j++) {
      fs.read((char*)m_inv_metric[i][j].data(), m_mesh.size() * sizeof(Scalar));
    }
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        fs.read(&m_connection_mask[i][j][k], sizeof(char));
        if (m_connection_mask[i][j][k] == 1)
          fs.read((char*)m_connection[i][j][k].data(), m_mesh.size() * sizeof(Scalar));
      }
    }
  }
  fs.close();
}

MultiArray<Scalar>&
Grid::scales(int n, Stagger_t stagger) {
  assert(n >= 0 && n < 3);
  return m_scales[n][stagger.to_ulong()];
}

const MultiArray<Scalar>&
Grid::scales(int n, Stagger_t stagger) const {
  assert(n >= 0 && n < 3);
  return m_scales[n][stagger.to_ulong()];
}

void
Grid::allocate_arrays() {
  m_det.resize(m_mesh.extent());
  m_alpha.resize(m_mesh.extent());
  for (int i = 0; i < 3; i++) {
    // m_cell_face_area[i].resize(m_mesh.extent());
    m_beta[i].resize(m_mesh.extent());
  }
}

Eigen::Matrix4d invert_matrix(const Eigen::Matrix4d& input) {
  // output = input.inverse();
  Eigen::Matrix4d result;
  result = input.inverse();
  return result;
}
std::array<std::array<double, 4>, 4>
Grid::invert_mat4(const std::array<std::array<double, 4>, 4>& mat) {
  std::array<std::array<double, 4>, 4> result;
  Eigen::Matrix4d n, m;
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      m(i, j) = mat[i][j];
  n = m.inverse();
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      result[i][j] = n(i, j);
  return result;
}
