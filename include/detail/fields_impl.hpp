#ifndef _FIELDS_IMPL_H_
#define _FIELDS_IMPL_H_

#include "fields.h"
#include "algorithms/interp_template.h"

namespace Aperture {

template <typename T>
template <typename Func>
void ScalarField<T>::initialize(const Func& f) {
  // This way scalar field is always defined in the center of the cell
  for (int k = 0; k < m_grid -> extent().depth(); ++k) {
    double x3 = m_grid -> mesh().pos(2, k, m_stagger[2]);
    for (int j = 0; j < m_grid -> extent().height(); ++j) {
      double x2 = m_grid -> mesh().pos(1, j, m_stagger[1]);
      for (int i = 0; i < m_grid -> extent().width(); ++i) {
        double x1 = m_grid -> mesh().pos(0, i, m_stagger[0]);
        m_array(i, j, k) = f(x1, x2, x3);
      }
    }
  }
}


template <typename T>
template <typename Func>
void VectorField<T>::initialize(int component, const Func& f) {
  // This way vector field is always defined in the center of the cell
  // face, staggered in the direction of the component
  for (int k = 0; k < m_grid -> extent().depth(); ++k) {
    double x3 = m_grid -> mesh().pos(2, k, m_stagger[component][2]);
    for (int j = 0; j < m_grid -> extent().height(); ++j) {
      double x2 = m_grid -> mesh().pos(1, j, m_stagger[component][1]);
      for (int i = 0; i < m_grid -> extent().width(); ++i) {
        double x1 = m_grid -> mesh().pos(0, i, m_stagger[component][0]);
        m_array[component](i, j, k) = f(x1, x2, x3);
      }
    }
  }
}

template <typename T>
template <typename Double>
Double
VectorField<T>::interpolate(int n, int cell, const Vec3<Double>& x) const {
  Double result = 0.0;
  auto mesh = m_grid -> mesh();

  Vec3<int> c = mesh.get_cell_3d(cell);
  Vec3<int> lower = c - Vec3<int>(1, 1, 1);
  Vec3<int> upper = c + Vec3<int>(1, 1, 1);
  if (m_grid -> dim() < 3) {
    lower[2] = upper[2] = c[2];
  }
  for (int k = lower[2]; k <= upper[2]; k++) {
    for (int j = lower[1]; j <= upper[1]; j++) {
      for (int i = lower[0]; i <= upper[0]; i++) {
        result += m_array[n](i, j, k) * interp_cell(x[0], i - c[0], mesh.delta[0], m_stagger[n][0])
                  * interp_cell(x[1], j - c[1], mesh.delta[1], m_stagger[n][1])
                  * (m_grid -> dim() < 3 ? 1.0 : interp_cell(x[2], k - c[2], mesh.delta[2], m_stagger[n][2]));
      }
    }
  }
  return result;
}

}

#endif  // _FIELDS_IMPL_H_
