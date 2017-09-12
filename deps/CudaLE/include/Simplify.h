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

// Generic single partial derivative class
template <int Argument, typename Expr>
struct Derivative;


template <typename Expr>
struct Simplified {
  typedef Expr arg_type;
  typedef Expr result_type;
  result_type result;

  HOST_DEVICE Simplified(arg_type expr) : result(expr) {}

  template <typename Data>
  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
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
typename std::enable_if<std::is_same<typename Simplified<Expr>::result_type, Expr>::value, Expr>::type
simplify(const Expr& expr) {
  // std::cout << "Cannot simplify further" << std::endl;
  // expr.print(); printf("\n");
  return expr;
}

template <typename Expr>
HD_INLINE
typename std::enable_if<std::is_same<typename Simplified<Expr>::result_type, Expr>::value == false,
                        typename Simplified<Expr>::result_type>::type
simplify(const Expr& expr) {
  // std::cout << "Can simplify further" << std::endl;
  // expr.print(); printf("\n");
  return Simplified<Expr>(expr).result;
}

// ////////////////////////////////////////////////////////////////////////////////
// ///  Operators
// ////////////////////////////////////////////////////////////////////////////////

template <typename Op, typename Left, typename Right>
struct Simplified<BinaryOp<Op, Left, Right> > {
  typedef BinaryOp<Op, Left, Right> arg_type;
  typedef BinaryOp<Op, typename Simplified<Left>::result_type,
                   typename Simplified<Right>::result_type > result_type;
  result_type result;

  HOST_DEVICE Simplified(arg_type expr) : result(simplify(expr.left), simplify(expr.right)) {}

  template <typename Data>
  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return result(x1, x2, x3, x4);
  }
};

template <typename Op, typename Arg>
struct Simplified<UnaryOp<Op, Arg> > {
  typedef UnaryOp<Op, Arg> arg_type;
  typedef UnaryOp<Op, typename Simplified<Arg>::result_type> result_type;
  result_type result;

  HOST_DEVICE Simplified(arg_type expr) : result(simplify(expr.arg)) {}

  template <typename Data>
  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return result(x1, x2, x3, x4);
  }
};

template <int Argument, typename Expr>
struct Simplified<Derivative<Argument, Expr>> {
  typedef Derivative<Argument, Expr> arg_type;
  typedef typename Simplified<typename Derivative<Argument, Expr>::result_type>::result_type result_type;
  result_type result;

  HOST_DEVICE Simplified(arg_type expr) : result(simplify(expr.derivative)) {}

  template <typename Data>
  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return result(x1, x2, x3, x4);
  }
};


////////////////////////////////////////////////////////////////////////////////
///  Explicit templates
////////////////////////////////////////////////////////////////////////////////

// Simplifying left + zero = left
template <typename Left>
struct Simplified<BinaryOp<Plus, Left, ZeroOp> > {
  typedef BinaryOp<Plus, Left, ZeroOp> arg_type;
  typedef typename Simplified<Left>::result_type result_type;
  result_type result;

  HOST_DEVICE Simplified(arg_type expr) : result(simplify(expr.left)) {}

  template <typename Data>
  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return result(x1, x2, x3, x4);
  }
};

// Simplifying zero + right = right
template <typename Right>
struct Simplified<BinaryOp<Plus, ZeroOp, Right> > {
  typedef BinaryOp<Plus, ZeroOp, Right> arg_type;
  typedef typename Simplified<Right>::result_type result_type;
  result_type result;

  HOST_DEVICE Simplified(arg_type expr) : result(simplify(expr.right)) {}

  template <typename Data>
  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return result(x1, x2, x3, x4);
  }
};

// Special case to resolve ambiguity
template <>
struct Simplified<BinaryOp<Plus, ZeroOp, ZeroOp> > {
  typedef BinaryOp<Plus, ZeroOp, ZeroOp> arg_type;
  typedef ZeroOp result_type;
  result_type result;

  HOST_DEVICE Simplified(arg_type expr) : result(ZeroOp{}) {}

  template <typename Data>
  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return result(x1, x2, x3, x4);
  }
};

// Simplifying left - zero
template <typename Left>
struct Simplified<BinaryOp<Minus, Left, ZeroOp> > {
  typedef BinaryOp<Minus, Left, ZeroOp> arg_type;
  typedef Left result_type;
  result_type result;

  HOST_DEVICE Simplified(arg_type expr) : result(expr.left) {}

  template <typename Data>
  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return result(x1, x2, x3, x4);
  }
};

// // Simplifying zero - right
// template <typename Right>
// struct Simplified<BinaryOp<Minus, ZeroOp, Right> > {
//   typedef BinaryOp<Minus, ZeroOp, Right> arg_type;
//   typedef typename Simplified<Right>::result_type result_type;
//   result_type result;

//   HOST_DEVICE Simplified(arg_type expr) : result(simplify(expr.right)) {}

//   template <typename Data>
//   HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
//     return result(x1, x2, x3, x4);
//   }
// };

// // Special case to resolve ambiguity
// template <>
// struct Simplified<BinaryOp<Minus, ZeroOp, ZeroOp> > {
//   typedef BinaryOp<Minus, ZeroOp, ZeroOp> arg_type;
//   typedef ZeroOp result_type;
//   result_type result;

//   HOST_DEVICE Simplified(arg_type expr) : result(ZeroOp{}) {}

//   template <typename Data>
//   HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
//     return result(x1, x2, x3, x4);
//   }
// };

// Simplify left * zero = zero
template <typename Left>
struct Simplified<BinaryOp<Multiply, Left, ZeroOp> > {
  typedef BinaryOp<Multiply, Left, ZeroOp> arg_type;
  typedef ZeroOp result_type;
  result_type result;

  HOST_DEVICE Simplified(arg_type expr) : result(ZeroOp{}) {}

  template <typename Data>
  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return result(x1, x2, x3, x4);
  }
};

// Simplify zero * right = zero
template <typename Right>
struct Simplified<BinaryOp<Multiply, ZeroOp, Right> > {
  typedef BinaryOp<Multiply, ZeroOp, Right> arg_type;
  typedef ZeroOp result_type;
  result_type result;

  HOST_DEVICE Simplified(arg_type expr) : result(ZeroOp{}) {}

  template <typename Data>
  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return result(x1, x2, x3, x4);
  }
};

// Simplify left * one = left
template <typename Left>
struct Simplified<BinaryOp<Multiply, Left, OneOp> > {
  typedef BinaryOp<Multiply, Left, OneOp> arg_type;
  typedef Left result_type;
  result_type result;

  HOST_DEVICE Simplified(arg_type expr) : result(expr.left) {}

  template <typename Data>
  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return result(x1, x2, x3, x4);
  }
};

// Simplify one * right = right
template <typename Right>
struct Simplified<BinaryOp<Multiply, OneOp, Right> > {
  typedef BinaryOp<Multiply, OneOp, Right> arg_type;
  typedef Right result_type;
  result_type result;

  HOST_DEVICE Simplified(arg_type expr) : result(expr.right) {}

  template <typename Data>
  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return result(x1, x2, x3, x4);
  }
};

// Ambiguous case of one * one = one
template <>
struct Simplified<BinaryOp<Multiply, OneOp, OneOp> > {
  typedef BinaryOp<Multiply, OneOp, OneOp> arg_type;
  typedef OneOp result_type;
  result_type result;

  HOST_DEVICE Simplified(arg_type expr) : result(OneOp{}) {}

  template <typename Data>
  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return result(x1, x2, x3, x4);
  }
};

// Simplify zero / right = zero
template <typename Right>
struct Simplified<BinaryOp<Divide, ZeroOp, Right> > {
  typedef BinaryOp<Divide, ZeroOp, Right> arg_type;
  typedef ZeroOp result_type;
  result_type result;

  HOST_DEVICE Simplified(arg_type expr) : result(ZeroOp{}) {}

  template <typename Data>
  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return result(x1, x2, x3, x4);
  }
};

// Simplify left / one = left
template <typename Left>
struct Simplified<BinaryOp<Divide, Left, OneOp> >
{
  typedef BinaryOp<Divide, Left, OneOp> arg_type;
  typedef Left result_type;
  result_type result;

  HOST_DEVICE Simplified(arg_type expr) : result(expr.left) {}

  template <typename Data>
  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return result(x1, x2, x3, x4);
  }
};

}

#endif  // _SIMPLIFY_H_
