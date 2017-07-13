/////////////////////////////////////////////////////////////////////////////////////////
///
///           \file  Simplify.h
///
/// __Description__:     Provides the simplify function to reduce the terms
///
/// __Last modified__:   <>\n
/// __Version__:         1.0\n
/// __Author__:          Alex Chen, fizban007@gmail.com\n
/// __Organization__:    Columbia University
///
/////////////////////////////////////////////////////////////////////////////////////////

#ifndef _SIMPLIFY_H_
#define _SIMPLIFY_H_

#include "Operators.h"
#include "Functions.h"

namespace CudaLE {

template <typename Expr>
struct Simplified {
  typedef Expr arg_type;
  typedef Expr result_type;
  result_type result;

  HOST_DEVICE Simplified(arg_type expr) : result(expr) {}

  HD_INLINE double operator() (double x1, double x2 = 0.0, double x3 = 0.0, double x4 = 0.0) {
    return result(x1, x2, x3, x4);
  }

  HD_INLINE void print() const {
    result.print();
  }
};

////////////////////////////////////////////////////////////////////////////////
///  Functions
////////////////////////////////////////////////////////////////////////////////

template <typename Expr>
HD_INLINE
typename Simplified<Expr>::result_type
simplify(const Expr& expr) {
  return Simplified<Expr>(expr).result;
}

////////////////////////////////////////////////////////////////////////////////
///  Explicit templates
////////////////////////////////////////////////////////////////////////////////

// Simplifying left + zero
template <typename Left>
struct Simplified<BinaryOp<Plus, Left, ZeroOp> > {
  typedef BinaryOp<Plus, Left, ZeroOp> arg_type;
  typedef typename Simplified<Left>::result_type result_type;
  result_type result;

  HOST_DEVICE Simplified(arg_type expr) : result(simplify(expr.left)) {}

  HD_INLINE double operator() (double x1, double x2 = 0.0, double x3 = 0.0, double x4 = 0.0) {
    return result(x1, x2, x3, x4);
  }
};

// Simplifying zero + right
template <typename Right>
struct Simplified<BinaryOp<Plus, ZeroOp, Right> > {
  typedef BinaryOp<Plus, ZeroOp, Right> arg_type;
  typedef typename Simplified<Right>::result_type result_type;
  result_type result;

  HOST_DEVICE Simplified(arg_type expr) : result(simplify(expr.right)) {}

  HD_INLINE double operator() (double x1, double x2 = 0.0, double x3 = 0.0, double x4 = 0.0) {
    return result(x1, x2, x3, x4);
  }
};

// Special case
template <>
struct Simplified<BinaryOp<Plus, ZeroOp, ZeroOp> > {
  typedef BinaryOp<Plus, ZeroOp, ZeroOp> arg_type;
  typedef ZeroOp result_type;
  result_type result;

  HOST_DEVICE Simplified(arg_type expr) : result(ZeroOp{}) {}

  HD_INLINE double operator() (double x1, double x2 = 0.0, double x3 = 0.0, double x4 = 0.0) {
    return result(x1, x2, x3, x4);
  }
};

// Simplifying left - zero
template <typename Left>
struct Simplified<BinaryOp<Minus, Left, ZeroOp> > {
  typedef BinaryOp<Minus, Left, ZeroOp> arg_type;
  typedef typename Simplified<Left>::result_type result_type;
  result_type result;

  HOST_DEVICE Simplified(arg_type expr) : result(simplify(expr.left)) {}

  HD_INLINE double operator() (double x1, double x2 = 0.0, double x3 = 0.0, double x4 = 0.0) {
    return result(x1, x2, x3, x4);
  }
};

// Simplifying zero - right
template <typename Right>
struct Simplified<BinaryOp<Minus, ZeroOp, Right> > {
  typedef BinaryOp<Minus, ZeroOp, Right> arg_type;
  typedef typename Simplified<Right>::result_type result_type;
  result_type result;

  HOST_DEVICE Simplified(arg_type expr) : result(simplify(expr.right)) {}

  HD_INLINE double operator() (double x1, double x2 = 0.0, double x3 = 0.0, double x4 = 0.0) {
    return result(x1, x2, x3, x4);
  }
};

// Special case
template <>
struct Simplified<BinaryOp<Minus, ZeroOp, ZeroOp> > {
  typedef BinaryOp<Minus, ZeroOp, ZeroOp> arg_type;
  typedef ZeroOp result_type;
  result_type result;

  HOST_DEVICE Simplified(arg_type expr) : result(ZeroOp{}) {}

  HD_INLINE double operator() (double x1, double x2 = 0.0, double x3 = 0.0, double x4 = 0.0) {
    return result(x1, x2, x3, x4);
  }
};

// Simplify left * zero
template <typename Left>
struct Simplified<BinaryOp<Multiply, Left, ZeroOp> > {
  typedef BinaryOp<Multiply, Left, ZeroOp> arg_type;
  typedef ZeroOp result_type;
  result_type result;

  HOST_DEVICE Simplified(arg_type expr) : result(ZeroOp{}) {}

  HD_INLINE double operator() (double x1, double x2 = 0.0, double x3 = 0.0, double x4 = 0.0) {
    return result(x1, x2, x3, x4);
  }
};

// Simplify left * zero
template <typename Right>
struct Simplified<BinaryOp<Multiply, ZeroOp, Right> > {
  typedef BinaryOp<Multiply, ZeroOp, Right> arg_type;
  typedef ZeroOp result_type;
  result_type result;

  HOST_DEVICE Simplified(arg_type expr) : result(ZeroOp{}) {}

  HD_INLINE double operator() (double x1, double x2 = 0.0, double x3 = 0.0, double x4 = 0.0) {
    return result(x1, x2, x3, x4);
  }
};

////////////////////////////////////////////////////////////////////////////////
///  Operators
////////////////////////////////////////////////////////////////////////////////

template <typename Left, typename Right>
struct Simplified<BinaryOp<Multiply, Left, Right> > {
  typedef BinaryOp<Multiply, Left, Right> arg_type;
  typedef BinaryOp<Multiply, typename Simplified<Left>::result_type,
                   typename Simplified<Right>::result_type > result_type;
  result_type result;

  HOST_DEVICE Simplified(arg_type expr) : result(simplify(expr.left), simplify(expr.right)) {}

  HD_INLINE double operator() (double x1, double x2 = 0.0, double x3 = 0.0, double x4 = 0.0) {
    return result(x1, x2, x3, x4);
  }
};

template <typename Left, typename Right>
struct Simplified<BinaryOp<Plus, Left, Right> > {
  typedef BinaryOp<Plus, Left, Right> arg_type;
  typedef BinaryOp<Plus, typename Simplified<Left>::result_type,
                   typename Simplified<Right>::result_type > result_type;
  result_type result;

  HOST_DEVICE Simplified(arg_type expr) : result(simplify(expr.left), simplify(expr.right)) {}

  HD_INLINE double operator() (double x1, double x2 = 0.0, double x3 = 0.0, double x4 = 0.0) {
    return result(x1, x2, x3, x4);
  }
};

template <typename Left, typename Right>
struct Simplified<BinaryOp<Minus, Left, Right> > {
  typedef BinaryOp<Minus, Left, Right> arg_type;
  typedef BinaryOp<Minus, typename Simplified<Left>::result_type,
                   typename Simplified<Right>::result_type > result_type;
  result_type result;

  HOST_DEVICE Simplified(arg_type expr) : result(simplify(expr.left), simplify(expr.right)) {}

  HD_INLINE double operator() (double x1, double x2 = 0.0, double x3 = 0.0, double x4 = 0.0) {
    return result(x1, x2, x3, x4);
  }
};

}

#endif  // _SIMPLIFY_H_
