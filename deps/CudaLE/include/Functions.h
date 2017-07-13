#ifndef  _FUNCTIONS_H_
#define  _FUNCTIONS_H_

#include <math.h>
#include "cudaControl.h"
#include "Operators.h"
#include "helper.h"

namespace CudaLE {

struct Sin
{
  HD_INLINE static double apply(double x) {
    return sin(x);
  }

  HD_INLINE static void print() {
    helper::print("sin");
  }
};

template <typename Arg>
HD_INLINE UnaryOp<Sin, Arg> sin(Arg arg) {
  return UnaryOp<Sin, Arg>(arg);
}

struct Cos
{
  HD_INLINE static double apply(double x) {
    return cos(x);
  }

  HD_INLINE static void print() {
    helper::print("cos");
  }
};

template <typename Arg>
HD_INLINE UnaryOp<Cos, Arg> cos(Arg arg) {
  return UnaryOp<Cos, Arg>(arg);
}

struct Exp
{
  HD_INLINE static double apply(double x) {
    return exp(x);
  }

  HD_INLINE static void print() {
    helper::print("exp");
  }
};

template <typename Arg>
HD_INLINE UnaryOp<Exp, Arg> exp(Arg arg) {
  return UnaryOp<Exp, Arg>(arg);
}

struct Log
{
  HD_INLINE static double apply(double x) {
    return log(x);
  }

  HD_INLINE static void print() {
    helper::print("log");
  }
};

template <typename Arg>
HD_INLINE UnaryOp<Log, Arg> log(Arg arg) {
  return UnaryOp<Log, Arg>(arg);
}

struct Plus
{
  HD_INLINE static double apply(double a, double b) {
    return a + b;
  }

  HD_INLINE static void print() {
    helper::print(" + ");
  }
};

template <typename Left, typename Right>
HD_INLINE BinaryOp<Plus, Left, Right> operator+ (Left lhs, Right rhs) {
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
  HD_INLINE static double apply(double a, double b) {
    return a - b;
  }

  HD_INLINE static void print() {
    helper::print(" - ");
  }
};

template <typename Left, typename Right>
HD_INLINE BinaryOp<Minus, Left, Right> operator- (Left lhs, Right rhs) {
  return BinaryOp<Minus, Left, Right>(lhs, rhs);
}

struct Multiply
{
  HD_INLINE static double apply(double a, double b) {
    return a * b;
  }

  HD_INLINE static void print() {
    helper::print(" * ");
  }
};

template <typename Left, typename Right>
HD_INLINE BinaryOp<Multiply, Left, Right> operator* (Left lhs, Right rhs) {
  return BinaryOp<Multiply, Left, Right>(lhs, rhs);
}

// template <typename Left>
// BinaryOp<Multiply, Left, ConstOp> operator* (Left lhs, double rhs) {
//     return BinaryOp<Multiply, Left, ConstOp>(lhs, ConstOp(rhs));
// }

// template <typename Right>
// BinaryOp<Multiply, ConstOp, Right> operator* (double lhs, Right rhs) {
//     return BinaryOp<Multiply, ConstOp, Right>(ConstOp(lhs), rhs);
// }

struct Divide
{
  HD_INLINE static double apply(double a, double b) {
    return a / b;
  }

  HD_INLINE static void print() {
    helper::print(" / ");
  }
};

template <typename Left, typename Right>
HD_INLINE BinaryOp<Divide, Left, Right> operator/ (Left lhs, Right rhs) {
  return BinaryOp<Divide, Left, Right>(lhs, rhs);
}

// template <typename Left>
// BinaryOp<Divide, Left, ConstOp> operator/ (Left lhs, double rhs) {
//     return BinaryOp<Divide, Left, ConstOp>(lhs, ConstOp(rhs));
// }

// template <typename Right>
// BinaryOp<Divide, ConstOp, Right> operator/ (double lhs, Right rhs) {
//     return BinaryOp<Divide, ConstOp, Right>(ConstOp(lhs), rhs);
// }

template<int power>
struct Pow
{
  HD_INLINE static double apply(double a) {
    return pow(a, power);
  }

  HD_INLINE static void print() {
    helper::print("pow<");
    helper::print(power);
    helper::print(">");
  }
};

template <int power, typename Arg>
HD_INLINE UnaryOp<Pow<power>, Arg> pow(Arg arg) {
  return UnaryOp<Pow<power>, Arg>(arg);
}

struct Sqrt
{
  HD_INLINE static double apply(double a) {
    return sqrt(a);
  }

  HD_INLINE static void print() {
    helper::print("sqrt");
  }
};

template <typename Arg>
HD_INLINE UnaryOp<Sqrt, Arg> sqrt(Arg arg) {
  return UnaryOp<Sqrt, Arg>(arg);
}

struct Square
{
  HD_INLINE static double apply(double a) {
    return a * a;
  }

  HD_INLINE static void print() {
    helper::print("square");
  }
};

template <typename Arg>
HD_INLINE UnaryOp<Square, Arg> square(Arg arg) {
  return UnaryOp<Square, Arg>(arg);
}

}



#endif   // ----- #ifndef _FUNCTIONS_H_  ----- 
