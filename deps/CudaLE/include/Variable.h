/////////////////////////////////////////////////////////////////////////////////////////
///
///           \file  Variable.h
///
/// __Description__:     Declares the variables used for the expression template
///
/// __Last modified__:   <2014-07-08 19:14:47 alex>\n
/// __Version__:         1.0\n
/// __Author__:          Alex Chen, fizban007@gmail.com\n
/// __Organization__:    Columbia University
///
/////////////////////////////////////////////////////////////////////////////////////////

#ifndef  _VARIABLE_H_
#define  _VARIABLE_H_

#include <type_traits>
#include "cudaControl.h"
#include "helper.h"

namespace CudaLE {

template <int Argument, typename Data = double>
struct Var
{
  typedef Var<Argument, Data> type;

  HOST_DEVICE Var() {}
  HD_INLINE Data operator() (Data x1, Data x2 = 0.0, Data x3 = 0.0, Data x4 = 0.0) const {
    if (1 == Argument)
      return x1;
    else if (2 == Argument)
      return x2;
    else if (3 == Argument)
      return x3;
    else if (4 == Argument)
      return x4;
    else
      return Data();
  }

  template <typename T>
  HD_INLINE
  bool
  operator== (const T& obj) const {
    return false;
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value, bool>::type
  operator== (const T& obj) const {
    return true;
  }

  template <typename T>
  HD_INLINE bool operator!= (const T& obj) const {
    return (!operator==(obj));
  }

  HD_INLINE void print() const {
    if (1 == Argument)
      helper::print("x1");
    else if (2 == Argument)
      helper::print("x2");
    else if (3 == Argument)
      helper::print("x3");
    else if (4 == Argument)
      helper::print("x4");
    else
      helper::print("");
  }
};

namespace placeholders {

static Var<1, double> _1;
static Var<2, double> _2;
static Var<3, double> _3;
static Var<4, double> _4;

namespace cartesian {
static Var<1, double> _x;
static Var<2, double> _y;
static Var<3, double> _z;
}

namespace spherical {
static Var<1, double> _r;
static Var<2, double> _theta;
static Var<3, double> _phi;
}

}

}

// #define _1 Var<1, double>()
// #define _2 Var<2, double>()
// #define _3 Var<3, double>()

#endif   // ----- #ifndef _VARIABLE_H_  ----- 
