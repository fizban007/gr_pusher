#include "metrics-analytic.h"
#include "fadiff.h"

typedef fadbad::F<double, 6> var;

namespace Aperture {

template struct MetricBase<var>;
template struct Schwarzschild<var>;
template struct Boyer_Lindquist<var>;
template struct Cartesian_KS<var>;
template struct Kerr_Schild<var>;

}
