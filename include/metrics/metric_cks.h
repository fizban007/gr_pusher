#ifndef _METRIC_CKS_H_
#define _METRIC_CKS_H_

#include "CudaLE.h"
#include "vec3.h"
#include "ptc.h"


namespace Aperture {

template <typename Data>
struct Cartesian_KS {
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

  DEFINE_FUNCTOR(g00, -1.0 + f);
  DEFINE_FUNCTOR(g01, f* l1);
  DEFINE_FUNCTOR(g02, f* l2);
  DEFINE_FUNCTOR(g03, f* l3);
  DEFINE_FUNCTOR(g10, f* l1);
  DEFINE_FUNCTOR(g11, 1.0 + f * square(l1));
  DEFINE_FUNCTOR(g12, f* l1* l2);
  DEFINE_FUNCTOR(g13, f* l1* l3);
  DEFINE_FUNCTOR(g20, f* l2);
  DEFINE_FUNCTOR(g21, f* l1* l2);
  DEFINE_FUNCTOR(g22, 1.0 + f * square(l2));
  DEFINE_FUNCTOR(g23, f* l2* l3);
  DEFINE_FUNCTOR(g30, f* l3);
  DEFINE_FUNCTOR(g31, f* l1* l3);
  DEFINE_FUNCTOR(g32, f* l2* l3);
  DEFINE_FUNCTOR(g33, 1.0 + f * square(l3));

  DEFINE_FUNCTOR(inv_g00, -1.0 - f);
  DEFINE_FUNCTOR(inv_g01, f* l1);
  DEFINE_FUNCTOR(inv_g02, f* l2);
  DEFINE_FUNCTOR(inv_g03, f* l3);
  DEFINE_FUNCTOR(inv_g10, f* l1);
  DEFINE_FUNCTOR(inv_g11, 1.0 - f * square(l1));
  DEFINE_FUNCTOR(inv_g12, -1.0 * f * l1 * l2);
  DEFINE_FUNCTOR(inv_g13, -1.0 * f * l1 * l3);
  DEFINE_FUNCTOR(inv_g20, f* l2);
  DEFINE_FUNCTOR(inv_g21, -1.0 * f * l1 * l2);
  DEFINE_FUNCTOR(inv_g22, 1.0 - f * square(l2));
  DEFINE_FUNCTOR(inv_g23, -1.0 * f * l2 * l3);
  DEFINE_FUNCTOR(inv_g30, f* l3);
  DEFINE_FUNCTOR(inv_g31, -1.0 * f * l1 * l3);
  DEFINE_FUNCTOR(inv_g32, -1.0 * f * l2 * l3);
  DEFINE_FUNCTOR(inv_g33, 1.0 - f * square(l3));
};

template <typename Data>
Vec3<Data>
mid_point(const Vec3<Data>& x, const Vec3<double>& x0);

template <typename Data>
Data
quadratic_solve(const Data& a, const Data& b, const Data& c);

template <typename Data>
Data
Gamma(const Vec3<Data>& x, const Vec3<Data>& u,
      const Cartesian_KS<Data>& metric);

template <typename Data>
Data
u0_energy(const Vec3<Data>& x, const Vec3<Data>& u,
          const Cartesian_KS<Data>& metric);

// template <typename Data>
// Data
// connection(int i, int a, int b, const Vec3<Data>& x,
//            const Cartesian_KS<Data>& metric) {
//   Data result;
//   if (i == 1) {

//   }
//   return result;
// }

template <typename Data>
Data
Fu(int n, const Vec3<double>& u0, const Vec3<double>& x0, const Vec3<Data>& u,
   const Vec3<Data>& x, const Cartesian_KS<Data>& metric, double dt);

template <typename Data>
Data
Fx(int n, const Vec3<double>& u0, const Vec3<double>& x0, const Vec3<Data>& u,
   const Vec3<Data>& x, const Cartesian_KS<Data>& metric, double dt);

template <typename Metric>
int
iterate_newton(Particle& p, const Metric& metric, double dt);



}

#endif  // _METRIC_CKS_H_
