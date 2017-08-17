#ifndef  _FUNCTIONS_H_
#define  _FUNCTIONS_H_

#include <math.h>
#include "cudaControl.h"
#include "Operators.h"
#include "helper.h"

namespace CudaLE {

struct Sin
{
  template <typename Data>
  HD_INLINE static auto apply(const Data& x) {
    return sin(x);
  }

  HD_INLINE static void print() {
    helper::print("sin");
  }
};

template <typename Arg>
HD_INLINE UnaryOp<Sin, Arg> sin(const Arg& arg) {
  return UnaryOp<Sin, Arg>(arg);
}

struct Cos
{
  template <typename Data>
  HD_INLINE static auto apply(const Data& x) {
    return cos(x);
  }

  HD_INLINE static void print() {
    helper::print("cos");
  }
};

template <typename Arg>
HD_INLINE UnaryOp<Cos, Arg> cos(const Arg& arg) {
  return UnaryOp<Cos, Arg>(arg);
}

struct Exp
{
  template <typename Data>
  HD_INLINE static auto apply(const Data& x) {
    return exp(x);
  }

  HD_INLINE static void print() {
    helper::print("exp");
  }
};

template <typename Arg>
HD_INLINE UnaryOp<Exp, Arg> exp(const Arg& arg) {
  return UnaryOp<Exp, Arg>(arg);
}

struct Log
{
  template <typename Data>
  HD_INLINE static auto apply(const Data& x) {
    return log(x);
  }

  HD_INLINE static void print() {
    helper::print("log");
  }
};

template <typename Arg>
HD_INLINE UnaryOp<Log, Arg> log(const Arg& arg) {
  return UnaryOp<Log, Arg>(arg);
}

struct Plus
{
  template <typename Left, typename Right>
  HD_INLINE static auto apply(const Left& a, const Right& b) {
    return a + b;
  }

  HD_INLINE static void print() {
    helper::print(" + ");
  }
};

template <typename Left, typename Right>
HD_INLINE BinaryOp<Plus, Left, Right> operator+ (const Left& lhs, const Right& rhs) {
  return BinaryOp<Plus, Left, Right>(lhs, rhs);
}

// template <typename Left>
// BinaryOp<Plus, Left, ConstOp> operator+ (Left lhs, double rhs) {
//     return BinaryOp<Plus, Left, ConstOp>(lhs, ConstOp(rhs));
// }

// template <typename Right>
// BinaryOp<Plus, ConstOp, Right> operator+ (double lhs, Right rhs) {
//     return BinaryOp<Plus, ConstOp, Right>(ConstOp(lhs), rhs);
// }

struct Minus
{
  template <typename Left, typename Right>
  HD_INLINE static auto apply(const Left& a, const Right& b) {
    return a - b;
  }

  HD_INLINE static void print() {
    helper::print(" - ");
  }
};

template <typename Left, typename Right>
HD_INLINE BinaryOp<Minus, Left, Right> operator- (const Left& lhs, const Right& rhs) {
  return BinaryOp<Minus, Left, Right>(lhs, rhs);
}

struct Multiply
{
  template <typename Left, typename Right>
  HD_INLINE static auto apply(const Left& a, const Right& b) {
    return a * b;
  }

  HD_INLINE static void print() {
    helper::print(" * ");
  }
};

template <typename Left, typename Right>
HD_INLINE BinaryOp<Multiply, Left, Right> operator* (const Left& lhs, const Right& rhs) {
  return BinaryOp<Multiply, Left, Right>(lhs, rhs);
}

struct Divide
{
  template <typename Left, typename Right>
  HD_INLINE static auto apply(const Left& a, const Right& b) {
    return a / b;
  }

  HD_INLINE static void print() {
    helper::print(" / ");
  }
};

template <typename Left, typename Right>
HD_INLINE BinaryOp<Divide, Left, Right> operator/ (const Left& lhs, const Right& rhs) {
  return BinaryOp<Divide, Left, Right>(lhs, rhs);
}

template<int power>
struct Pow
{
  template <typename Data>
  HD_INLINE static auto apply(const Data& a) {
    return pow(a, power);
  }

  HD_INLINE static void print() {
    helper::print("pow<");
    helper::print(power);
    helper::print(">");
  }
};

template <int power, typename Arg>
HD_INLINE UnaryOp<Pow<power>, Arg> pow(const Arg& arg) {
  return UnaryOp<Pow<power>, Arg>(arg);
}

struct Sqrt
{
  template <typename Data>
  HD_INLINE static auto apply(const Data& a) {
    return sqrt(a);
  }

  HD_INLINE static void print() {
    helper::print("sqrt");
  }
};

template <typename Arg>
HD_INLINE UnaryOp<Sqrt, Arg> sqrt(const Arg& arg) {
  return UnaryOp<Sqrt, Arg>(arg);
}

struct Square
{
  template <typename Data>
  HD_INLINE static auto apply(const Data& a) {
    return a * a;
  }

  HD_INLINE static void print() {
    helper::print("square");
  }
};

template <typename Arg>
HD_INLINE UnaryOp<Square, Arg> square(const Arg& arg) {
  return UnaryOp<Square, Arg>(arg);
}

}



#endif   // ----- #ifndef _FUNCTIONS_H_  -----
