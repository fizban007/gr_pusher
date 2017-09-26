#ifndef _METRICS_ANALYTIC_H_
#define _METRICS_ANALYTIC_H_

#include "CudaLE.h"

namespace Aperture {

using CudaLE::ZeroOp;
using CudaLE::OneOp;

template <typename Data>
struct MetricBase {
  DEFINE_FUNCTOR(g00, CudaLE::ConstOp(-1.0));

  DEFINE_FUNCTOR(g01, ZeroOp{});
  DEFINE_FUNCTOR(g02, ZeroOp{});
  DEFINE_FUNCTOR(g03, ZeroOp{});
  DEFINE_FUNCTOR(g10, ZeroOp{});
  DEFINE_FUNCTOR(g20, ZeroOp{});
  DEFINE_FUNCTOR(g30, ZeroOp{});

  DEFINE_FUNCTOR(g11, CudaLE::ConstOp(1.0));
  DEFINE_FUNCTOR(g22, CudaLE::ConstOp(1.0));
  DEFINE_FUNCTOR(g33, CudaLE::ConstOp(1.0));

  DEFINE_FUNCTOR(g12, ZeroOp{});
  DEFINE_FUNCTOR(g21, ZeroOp{});
  DEFINE_FUNCTOR(g13, ZeroOp{});
  DEFINE_FUNCTOR(g31, ZeroOp{});
  DEFINE_FUNCTOR(g23, ZeroOp{});
  DEFINE_FUNCTOR(g32, ZeroOp{});

  DEFINE_FUNCTOR(alpha, CudaLE::ConstOp(1.0));
  DEFINE_FUNCTOR(a2, CudaLE::ConstOp(1.0));
  DEFINE_FUNCTOR(b1, ZeroOp{});
  DEFINE_FUNCTOR(b2, ZeroOp{});
  DEFINE_FUNCTOR(b3, ZeroOp{});

  DEFINE_FUNCTOR(inv_g00, CudaLE::ConstOp(-1.0));

  DEFINE_FUNCTOR(inv_g01, ZeroOp{});
  DEFINE_FUNCTOR(inv_g02, ZeroOp{});
  DEFINE_FUNCTOR(inv_g03, ZeroOp{});
  DEFINE_FUNCTOR(inv_g10, ZeroOp{});
  DEFINE_FUNCTOR(inv_g20, ZeroOp{});
  DEFINE_FUNCTOR(inv_g30, ZeroOp{});

  DEFINE_FUNCTOR(inv_g11, CudaLE::ConstOp(1.0));
  DEFINE_FUNCTOR(inv_g22, CudaLE::ConstOp(1.0));
  DEFINE_FUNCTOR(inv_g33, CudaLE::ConstOp(1.0));

  DEFINE_FUNCTOR(inv_g12, ZeroOp{});
  DEFINE_FUNCTOR(inv_g21, ZeroOp{});
  DEFINE_FUNCTOR(inv_g13, ZeroOp{});
  DEFINE_FUNCTOR(inv_g31, ZeroOp{});
  DEFINE_FUNCTOR(inv_g23, ZeroOp{});
  DEFINE_FUNCTOR(inv_g32, ZeroOp{});
};

template <typename Data>
struct Schwarzschild : public MetricBase<Data> {
  CudaLE::Var<1, Data> _r;
  CudaLE::Var<2, Data> _theta;
  CudaLE::Var<3, Data> _phi;

  double m_ = 1.0;

  Schwarzschild() {}
  Schwarzschild(double m) : m_(m) {}

  DEFINE_FUNCTOR(g00, (2.0 * m_ / _r - 1.0));
  DEFINE_FUNCTOR(g11, 1.0 / (1.0 - 2.0 * m_ / _r));
  DEFINE_FUNCTOR(g22, _r* _r);
  DEFINE_FUNCTOR(g33, _r* _r* square(sin(_theta)));

  DEFINE_FUNCTOR(alpha, sqrt(1.0 - 2.0 * m_ / _r));
  DEFINE_FUNCTOR(a2, 1.0 - 2.0 * m_ / _r);

  DEFINE_FUNCTOR(inv_g00, 1.0 / (2.0 * m_ / _r - 1.0));
  DEFINE_FUNCTOR(inv_g11, 1.0 - 2.0 * m_ / _r);
  DEFINE_FUNCTOR(inv_g22, 1.0 / square(_r));
  DEFINE_FUNCTOR(inv_g33, 1.0 / square(_r * sin(_theta)));
};

template <typename Data>
struct Boyer_Lindquist : public MetricBase<Data> {
  CudaLE::Var<1, Data> _r;
  CudaLE::Var<2, Data> _theta;
  CudaLE::Var<3, Data> _phi;

  double m_ = 1.0;
  double a_ = 0.0;

  Boyer_Lindquist() {}
  Boyer_Lindquist(double m, double a) : m_(m), a_(a) {}

  // DEFINE_FUNCTOR(g00, (rg / _r - 1.0));
  // DEFINE_FUNCTOR(g11, 1.0 / (1.0 - rg / _r));
  // DEFINE_FUNCTOR(g22, _r* _r);
  // DEFINE_FUNCTOR(g33, _r* _r* square(sin(_theta)));

  DEFINE_FUNCTOR(rho2, _r* _r + a_ * a_ * square(cos(_theta)));
  DEFINE_FUNCTOR(z, 2.0 * m_ * _r / rho2);
  DEFINE_FUNCTOR(delta, _r* _r + a_ * a_ - 2.0 * m_ * _r);
  DEFINE_FUNCTOR(sigma,
                 square(_r* _r + a_ * a_) - delta * square(a_ * sin(_theta)));

  DEFINE_FUNCTOR(g00, z - 1.0);
  DEFINE_FUNCTOR(g11, rho2 / delta);
  DEFINE_FUNCTOR(g22, rho2);
  DEFINE_FUNCTOR(g33, sigma* square(sin(_theta)) / rho2);
  DEFINE_FUNCTOR(g03, -a_ * z * square(sin(_theta)));
  // DEFINE_FUNCTOR(g03, ZeroOp{});
  DEFINE_FUNCTOR(g30, -a_ * z * square(sin(_theta)));
  // DEFINE_FUNCTOR(g30, ZeroOp{});

  DEFINE_FUNCTOR(det, -1.0 * square(sin(_theta) * rho2));
  DEFINE_FUNCTOR(b3, -a_ * 2.0 * m_ * _r / sigma);
  DEFINE_FUNCTOR(a2, rho2* delta / sigma);
  DEFINE_FUNCTOR(alpha, sqrt(rho2* delta / sigma));

  DEFINE_FUNCTOR(inv_g00, -1.0 / a2);
  DEFINE_FUNCTOR(inv_g11, 1.0 / g11);
  DEFINE_FUNCTOR(inv_g22, 1.0 / g22);
  DEFINE_FUNCTOR(inv_g33, (1.0 / g33) - (square(b3) / a2));
  // DEFINE_FUNCTOR(inv_g33, 1.0 / g33);
  DEFINE_FUNCTOR(inv_g30, b3 / a2);
  // DEFINE_FUNCTOR(inv_g30, ZeroOp{});
  DEFINE_FUNCTOR(inv_g03, b3 / a2);
  // DEFINE_FUNCTOR(inv_g03, ZeroOp{});
};

template <typename Data>
struct Cartesian_KS : public MetricBase<Data> {
  CudaLE::Var<1, Data> _x;
  CudaLE::Var<2, Data> _y;
  CudaLE::Var<3, Data> _z;

  double m_ = 1.0;
  double a_ = 0.0;

  Cartesian_KS() {}
  Cartesian_KS(double m, double a) : m_(m), a_(a) {}

  // DEFINE_FUNCTOR(g00, (rg / _r - 1.0));
  // DEFINE_FUNCTOR(g11, 1.0 / (1.0 - rg / _r));
  // DEFINE_FUNCTOR(g22, _r* _r);
  // DEFINE_FUNCTOR(g33, _r* _r* square(sin(_theta)));

  DEFINE_FUNCTOR(R2, square(_x) + square(_y) + square(_z));
  DEFINE_FUNCTOR(r,
                 sqrt(0.5 * (R2 - a_ * a_ + sqrt(square(R2 - a_ * a_) +
                                                 4.0 * a_ * a_ * square(_z)))));
  DEFINE_FUNCTOR(r2,
                 0.5 * (R2 - a_ * a_ + sqrt(square(R2 - a_ * a_) +
                                            4.0 * a_ * a_ * square(_z))));
  DEFINE_FUNCTOR(l1, (r * _x + a_ * _y) / (r2 + a_ * a_));
  DEFINE_FUNCTOR(l2, (r * _y - a_ * _x) / (r2 + a_ * a_));
  DEFINE_FUNCTOR(l3, _z / r);
  DEFINE_FUNCTOR(f, 2.0 * m_ * r * r2 / (r2 * r2 + a_ * a_ * square(_z)));

  DEFINE_FUNCTOR(alpha, sqrt(OneOp{} / (OneOp{} + f)));
  DEFINE_FUNCTOR(a2, OneOp{} / (OneOp{} + f));

  DEFINE_FUNCTOR(g00, -1.0 + f);
  DEFINE_FUNCTOR(g01, f* l1);
  DEFINE_FUNCTOR(g02, f* l2);
  DEFINE_FUNCTOR(g03, f* l3);
  DEFINE_FUNCTOR(g10, f* l1);
  DEFINE_FUNCTOR(g11, OneOp{} + f * square(l1));
  DEFINE_FUNCTOR(g12, f* l1* l2);
  DEFINE_FUNCTOR(g13, f* l1* l3);
  DEFINE_FUNCTOR(g20, f* l2);
  DEFINE_FUNCTOR(g21, f* l1* l2);
  DEFINE_FUNCTOR(g22, OneOp{} + f * square(l2));
  DEFINE_FUNCTOR(g23, f* l2* l3);
  DEFINE_FUNCTOR(g30, f* l3);
  DEFINE_FUNCTOR(g31, f* l1* l3);
  DEFINE_FUNCTOR(g32, f* l2* l3);
  DEFINE_FUNCTOR(g33, OneOp{} + f * square(l3));

  DEFINE_FUNCTOR(inv_g00, -1.0 - f);
  DEFINE_FUNCTOR(inv_g01, f* l1);
  DEFINE_FUNCTOR(inv_g02, f* l2);
  DEFINE_FUNCTOR(inv_g03, f* l3);
  DEFINE_FUNCTOR(inv_g10, f* l1);
  DEFINE_FUNCTOR(inv_g11, OneOp{} - f * square(l1));
  DEFINE_FUNCTOR(inv_g12, -1.0 * f * l1 * l2);
  DEFINE_FUNCTOR(inv_g13, -1.0 * f * l1 * l3);
  DEFINE_FUNCTOR(inv_g20, f* l2);
  DEFINE_FUNCTOR(inv_g21, -1.0 * f * l1 * l2);
  DEFINE_FUNCTOR(inv_g22, OneOp{} - f * square(l2));
  DEFINE_FUNCTOR(inv_g23, -1.0 * f * l2 * l3);
  DEFINE_FUNCTOR(inv_g30, f* l3);
  DEFINE_FUNCTOR(inv_g31, -1.0 * f * l1 * l3);
  DEFINE_FUNCTOR(inv_g32, -1.0 * f * l2 * l3);
  DEFINE_FUNCTOR(inv_g33, OneOp{} - f * square(l3));

  DEFINE_FUNCTOR(b1, a2 * inv_g01);
  DEFINE_FUNCTOR(b2, a2 * inv_g02);
  DEFINE_FUNCTOR(b3, a2 * inv_g03);
};

template <typename Data>
struct Kerr_Schild : public MetricBase<Data> {
  CudaLE::Var<1, Data> _r;
  CudaLE::Var<2, Data> _theta;
  CudaLE::Var<3, Data> _phi;

  double m_ = 1.0;
  double a_ = 0.0;

  Kerr_Schild() {}
  Kerr_Schild(double m, double a) : m_(m), a_(a) {}

  DEFINE_FUNCTOR(rho2, _r * _r + a_ * a_ * square(cos(_theta)));
  DEFINE_FUNCTOR(z, 2.0 * m_ * _r / rho2);
  DEFINE_FUNCTOR(delta, _r * _r + a_ * a_ - 2.0 * m_ * _r);

  DEFINE_FUNCTOR(g00, z - 1.0);
  DEFINE_FUNCTOR(g11, 1.0 + z);
  DEFINE_FUNCTOR(g22, rho2);
  DEFINE_FUNCTOR(g33, (square(_r * _r + a_ * a_) -
                       square(a_ * sin(_theta)) * delta) *
                 square(sin(_theta)) / g22);
  DEFINE_FUNCTOR(g13, -a_ * square(sin(_theta)) * g11);
  DEFINE_FUNCTOR(g31, g13);

  DEFINE_FUNCTOR(alpha, 1.0 / sqrt(g11));
  DEFINE_FUNCTOR(a2, 1.0 / (1.0 + z));
  DEFINE_FUNCTOR(b1, z / (1.0 + z));

  DEFINE_FUNCTOR(g01, z);
  DEFINE_FUNCTOR(g03, -1.0 * z * a_ * square(sin(_theta)));
  DEFINE_FUNCTOR(g10, z);
  DEFINE_FUNCTOR(g30, -1.0 * z * a_ * square(sin(_theta)));
};

}

#endif  // _METRICS_ANALYTIC_H_
