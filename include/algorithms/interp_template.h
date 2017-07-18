#ifndef _INTERP_TEMPLATE_H_
#define _INTERP_TEMPLATE_H_

// #include <algorithm>

namespace Aperture {

template <typename Double>
Double
interp(const Double& dx) {
  // return max(1.0 - sqrt(sqr(dx)), 0.0);
  if (dx < 1.0 && dx >= 0.0) return 1.0 - dx;
  else if (dx < 0.0 && dx >= -1.0) return 1.0 + dx;
  else return 0.0;
}

template <typename Double>
Double
interp_cell(const Double& pos, int cell_diff, double delta, int stagger = 0) {
  // cell_diff = target_cell - ptc_cell
  Double x =
      (double)cell_diff + (stagger == 0 ? 0.5 : 1.0) - pos / delta;
  return interp(x);
}

template <typename Double>
Double
interp2(const Double& dx) {
  // return max(1.0 - sqrt(sqr(dx)), 0.0);
  // if (dx < 1.0 && dx >= 0.0) return 1.0 - dx;
  // else if (dx > -1.0 && dx < 0.0) return 1.0 + dx;
  if (dx < 2.0 && dx >= 1.0) return (2.0 - dx) * (2.0 - dx) * (2.0 - dx) / 6.0;
  else if (dx < 1.0 && dx >= 0.0) return (4.0 - 6.0 * dx * dx + 3.0 * dx * dx * dx) / 6.0;
  else if (dx < 0.0 && dx >= -1.0) return (4.0 - 6.0 * dx * dx - 3.0 * dx * dx * dx) / 6.0;
  if (dx < -1.0 && dx >= -2.0) return (2.0 + dx) * (2.0 + dx) * (2.0 + dx) / 6.0;
  else return 0.0;
}

template <typename Double>
Double
interp2_cell(const Double& pos, int cell_diff, double delta, int stagger = 0) {
  // cell_diff = target_cell - ptc_cell
  Double x =
      (double)cell_diff + (stagger == 0 ? 0.5 : 1.0) - pos / delta;
  return interp2(x);
}

}

#endif  // _INTERP_TEMPLATE_H_
