#include <iostream>
#include <math.h>
#include "CudaLE.h"

using namespace CudaLE;

using namespace CudaLE::placeholders;

// #define MY_DECLARE_FUNCTOR_TYPE(a, b) typedef typeof(a) b
// #define DEFINE_FUNCTOR(NAME, FUNCTOR)      \
//     typedef typeof(FUNCTOR) type_ ## NAME; \
//     type_ ## NAME NAME

class FuncTest
{
 public:
  DEFINE_FUNCTOR(f, cos(3.0 * _1) + exp(3.0 * _2));
  DEFINE_FUNCTOR(g, _1 * _2 * _2);
  DEFINE_FUNCTOR(h, pow<4>(_1));
  

  FuncTest() {}
  virtual ~FuncTest() {}
}; // ----- end of class FuncTest -----


int main()
{
  FuncTest F;
  // typeof (_1) func;
  // typedef typeof(_1 + _2) functype;
  // functype func = _1 + _2;
  // std::cout << func(2.0, 3.0) << std::endl;
  std::cout << F.f(2.0, 3.0) << std::endl;
  std::cout << (2.0 * F.f * F.g)(2.0, 3.0) << std::endl;
  std::cout << (D<1>(F.f) / D<2>(F.f))(2.0, 3.0) << " " << -3.0 * sin(6.0) << " " << 3.0 << std::endl;
  std::cout << D<2>(D<1>(F.g))(2.0, 3.0) << std::endl;
  std::cout << D<1>(F.h)(3.0) << std::endl;
  std::cout << D<2>(F.h)(3.0) << std::endl;
  helper::println(F.f);
  helper::println(F.g);
  helper::println(F.h);
  // std::cout << (ConstOp(3.0) * f.get_f ())(2.0, 3.0) << std::endl;
  // std::cout << func(2.0, 3.0) << std::endl;
  // std::cout << f.get_functor2()(2.0, 3.0) << std::endl;
  return 0;
}










