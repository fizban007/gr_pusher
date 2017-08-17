#ifndef  _DERIVATIVE_H_
#define  _DERIVATIVE_H_

#include <math.h>
#include "cudaControl.h"
#include "Functions.h"
#include "Variable.h"
#include "Simplify.h"

namespace CudaLE {

////////////////////////////////////////////////////////////////////////////////
///  Helper struct
////////////////////////////////////////////////////////////////////////////////
template<int> struct int2type{};

////////////////////////////////////////////////////////////////////////////////
///  Forward declaration
////////////////////////////////////////////////////////////////////////////////

// Generic single partial derivative class
template <int Argument, typename Expr>
struct Derivative;

// // Higher order partial derivative with respect to 1 argument
// template <int Argument, int times, typename Expr>
// struct Derivative_n;

// Higher order partial derivative with respect to 2 arguments
// template <int Arg1, int n1, int Arg2, int n2, typename Expr>
// struct Derivative_2n;

template <int Argument, int Arg2, typename Expr>
HD_INLINE
typename Derivative<Argument, typename Derivative<Arg2, Expr>::result_type>::result_type
D(const Derivative<Arg2, Expr>& deriv);

template <int Argument, typename Op, typename Left, typename Right>
HD_INLINE
typename Derivative<Argument, BinaryOp<Op, Left, Right> >::result_type
D(const BinaryOp<Op, Left, Right>& expr);

template <int Argument,typename Op, typename Right>
HD_INLINE
typename Derivative<Argument, BinaryOp<Op, double, Right> >::result_type
D(const BinaryOp<Op, double, Right>& expr);

template <int Argument,typename Op, typename Left>
HD_INLINE
typename Derivative<Argument, BinaryOp<Op, Left, double> >::result_type
D(const BinaryOp<Op, Left, double>& expr);

template <int Argument,typename Op, typename Expr>
HD_INLINE
typename Derivative<Argument, UnaryOp<Op, Expr> >::result_type
D(const UnaryOp<Op, Expr>& expr);

template <int Argument, int n, typename Data>
HD_INLINE
typename Derivative<Argument, Var<n, Data> >::result_type
D(const Var<n, Data>& var);

template <int Argument>
HD_INLINE
typename Derivative<Argument, double>::result_type
D(const double& val);

// Higher order derivatives
// template <int Argument, int times, typename Expr>
// HD_INLINE
// typename Derivative_n<Argument, times, Expr>::result_type
// D_n(const Expr& expr);

// template <int Arg1, int n1, int Arg2, int n2, typename Expr>
// HD_INLINE
// typename Derivative_2n<Arg1, n1, Arg2, n2, Expr>::result_type
// D_nm(const Expr& expr);

////////////////////////////////////////////////////////////////////////////////
///  Template classes
////////////////////////////////////////////////////////////////////////////////

// Derivative of a derivative
template <int Argument, int Arg2, typename Expr>
struct Derivative<Argument, Derivative<Arg2, Expr> >
{
  typedef typename Derivative<Arg2, Expr>::result_type arg_type;
  typedef typename Simplified<typename Derivative<Argument, arg_type>::result_type>::result_type result_type;
  result_type derivative;

  HOST_DEVICE Derivative(arg_type expr) : derivative(Simplify(D<Argument>(expr))) {}

  template <typename Data>
  HD_INLINE auto operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return derivative(x1, x2, x3, x4);
  }
};

// Derivative sum rule
template <int Argument, typename Left, typename Right>
struct Derivative<Argument, BinaryOp<Plus, Left, Right> >
{
  // Left left;
  // Right right;
  typedef BinaryOp<Plus, Left, Right> arg_type;
  typedef typename Simplified<
    BinaryOp<Plus, typename Derivative<Argument, Left>::result_type
             , typename Derivative<Argument, Right>::result_type > >::result_type result_type;
  result_type derivative;

  // HOST_DEVICE Derivative(Left t1, Right t2) : left(t1), right(t2),  {}
  HOST_DEVICE Derivative(Left t1, Right t2) : derivative(simplify(D<Argument>(t1) + D<Argument>(t2))) {}
  HOST_DEVICE Derivative(arg_type expr) : derivative(simplify(D<Argument>(expr.left) + D<Argument>(expr.right))) {}

  template <typename Data>
  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return derivative(x1, x2, x3, x4);
  }
};

// Derivative subtraction rule
template <int Argument, typename Left, typename Right>
struct Derivative<Argument, BinaryOp<Minus, Left, Right> >
{
  // Left left;
  // Right right;
  typedef BinaryOp<Minus, Left, Right> arg_type;
  typedef typename Simplified<
    BinaryOp<Minus, typename Derivative<Argument, Left>::result_type
             , typename Derivative<Argument, Right>::result_type > >::result_type result_type;
  result_type derivative;

  // HOST_DEVICE Derivative(Left t1, Right t2) : derivative(D<Argument>(t1), D<Argument>(t2))  {}
  HOST_DEVICE Derivative(Left t1, Right t2) : derivative(simplify(D<Argument>(t1) - D<Argument>(t2))) {}
  HOST_DEVICE Derivative(arg_type expr) : derivative(D<Argument>(expr.left), D<Argument>(expr.right)) {}
  HOST_DEVICE Derivative() {}

  template <typename Data>
  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return derivative(x1, x2, x3, x4);
  }
};

// Derivative multiplication rule
template <int Argument, typename Left, typename Right>
struct Derivative<Argument, BinaryOp<Multiply, Left, Right> >
{
  // Left left;
  // Right right;
  typedef BinaryOp<Multiply, Left, Right> arg_type;
  typedef typename Simplified<
    BinaryOp<Plus,
             typename Simplified<BinaryOp<Multiply, Left, typename Derivative<Argument, Right>::result_type> >::result_type,  typename Simplified<BinaryOp<Multiply, typename Derivative<Argument, Left>::result_type, Right> >::result_type > >::result_type result_type;
  result_type derivative;

  HOST_DEVICE Derivative(Left t1, Right t2) : derivative(simplify(simplify(t1 * D<Argument>(t2)) + simplify(D<Argument>(t1) * t2))) {}
  HOST_DEVICE Derivative(arg_type expr) :
      derivative(simplify(expr.left * D<Argument>(expr.right) + D<Argument>(expr.left) * expr.right)) {}
  HOST_DEVICE Derivative() {}

  template <typename Data>
  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return derivative(x1, x2, x3, x4);
  }
};

// Derivative division rule
template <int Argument, typename Left, typename Right>
struct Derivative<Argument, BinaryOp<Divide, Left, Right> >
{
  // Left left;
  // Right right;
  typedef BinaryOp<Divide, Left, Right> arg_type;
  typedef BinaryOp<Divide,
                   typename Simplified<BinaryOp<Minus,
                                    BinaryOp<Multiply, typename Derivative<Argument, Left>::result_type, Right>,
                                                BinaryOp<Multiply, Left, typename Derivative<Argument, Right>::result_type> > >::result_type,
                   BinaryOp<Multiply, Right, Right> > result_type;
  result_type derivative;

  HOST_DEVICE Derivative(Left t1, Right t2) :
      derivative(simplify(D<Argument>(t1) * t2 - t1 * D<Argument>(t2)), t2 * t2) {}
  HOST_DEVICE Derivative(arg_type expr) :
      derivative(D<Argument>(expr.left) * expr.right - expr.left * D<Argument>(expr.right), expr.right * expr.right) {}

  HOST_DEVICE Derivative() {}

  template <typename Data>
  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return derivative(x1, x2, x3, x4);
  }
};

// Derivative of sin function
template <int Argument, typename Arg>
struct Derivative<Argument, UnaryOp<Sin, Arg> >
{
  typedef UnaryOp<Sin, Arg> arg_type;
  typedef BinaryOp<Multiply, UnaryOp<Cos, Arg>,
                   typename Derivative<Argument, Arg>::result_type > result_type;
  result_type derivative;

  HOST_DEVICE Derivative(Arg arg) : derivative(cos(arg), D<Argument>(arg)) {}
  HOST_DEVICE Derivative(arg_type expr) : derivative(cos(expr.arg), D<Argument>(expr.arg)) {}

  template <typename Data>
  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return derivative(x1, x2, x3, x4);
  }
};

// Derivative of cos function
template <int Argument, typename Arg>
struct Derivative<Argument, UnaryOp<Cos, Arg> >
{
  typedef UnaryOp<Sin, Arg> arg_type;
  typedef BinaryOp<Multiply, BinaryOp<Multiply, double, UnaryOp<Sin, Arg> >,
                   typename Derivative<Argument, Arg>::result_type > result_type;
  result_type derivative;

  HOST_DEVICE Derivative(Arg arg) : derivative(-1.0 * sin(arg), D<Argument>(arg)) {}
  HOST_DEVICE Derivative(arg_type expr) : derivative(-1.0 * sin(expr.arg), D<Argument>(expr.arg)) {}

  template <typename Data>
  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return derivative(x1, x2, x3, x4);
  }
};

// Derivative of exponential function
template <int Argument, typename Arg>
struct Derivative<Argument, UnaryOp<Exp, Arg> >
{
  typedef UnaryOp<Exp, Arg> arg_type;
  typedef BinaryOp<Multiply, UnaryOp<Exp, Arg>,
                   typename Derivative<Argument, Arg>::result_type > result_type;
  result_type derivative;

  HOST_DEVICE Derivative(Arg arg) : derivative(exp(arg), D<Argument>(arg)) {}
  HOST_DEVICE Derivative(arg_type expr) : derivative(exp(expr.arg), D<Argument>(expr.arg)) {}

  template <typename Data>
  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return derivative(x1, x2, x3, x4);
  }
};

// Derivative of logarithm function
template <int Argument, typename Arg>
struct Derivative<Argument, UnaryOp<Log, Arg> >
{
  typedef UnaryOp<Log, Arg> arg_type;
  typedef BinaryOp<Divide, typename Derivative<Argument, Arg>::result_type, Arg> result_type;
  result_type derivative;

  HOST_DEVICE Derivative(Arg arg) : derivative(D<Argument>(arg), arg) {}
  HOST_DEVICE Derivative(arg_type expr) : derivative(D<Argument>(expr.arg), expr.arg) {}

  template <typename Data>
  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return derivative(x1, x2, x3, x4);
  }
};

// Derivative of an integer power function
template <int Argument, int n, typename Arg>
struct Derivative<Argument, UnaryOp<Pow<n>, Arg> >
{
  typedef UnaryOp<Pow<n>, Arg> arg_type;
  typedef BinaryOp<Multiply,
                   BinaryOp<Multiply, double, typename Derivative<Argument, Arg>::result_type>
                   , UnaryOp<Pow< n-1 >, Arg> >
      result_type;
  result_type derivative;

  HOST_DEVICE Derivative(Arg arg) : derivative((double)n * D<Argument>(arg), pow<n-1>(arg)) {}
  HOST_DEVICE Derivative(arg_type expr) : derivative((double)n * D<Argument>(expr.arg), pow<n-1>(expr.arg)) {}

  template <typename Data>
  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return derivative(x1, x2, x3, x4);
  }
};

// Derivative of square
template <int Argument, typename Arg>
struct Derivative<Argument, UnaryOp<Square, Arg> >
{
  typedef UnaryOp<Square, Arg> arg_type;
  typedef typename Simplified<BinaryOp<Multiply,
                                       typename Simplified<BinaryOp<Multiply, ConstOp, typename Derivative<Argument, Arg>::result_type> >::result_type
                                     , Arg> >::result_type
  result_type;
  result_type derivative;

  HOST_DEVICE Derivative(Arg arg) : derivative(simplify(simplify(ConstOp(2.0) * D<Argument>(arg)) * arg)) {}
  HOST_DEVICE Derivative(arg_type expr) : derivative(simplify(simplify(ConstOp(2.0) * D<Argument>(expr.arg)) * expr.arg)) {}

  template <typename Data>
  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return derivative(x1, x2, x3, x4);
  }
};

// Derivative of square root function
template <int Argument, typename Arg>
struct Derivative<Argument, UnaryOp<Sqrt, Arg> >
{
  typedef UnaryOp<Sqrt, Arg> arg_type;
  typedef BinaryOp<Divide,
                   typename Derivative<Argument, Arg>::result_type
                   , BinaryOp<Multiply, ConstOp, UnaryOp<Sqrt, Arg> > >
  result_type;
  result_type derivative;

  HOST_DEVICE Derivative(Arg arg) : derivative(D<Argument>(arg), ConstOp(2.0) * sqrt(arg)) {}
  HOST_DEVICE Derivative(arg_type expr) : derivative(D<Argument>(expr.arg), ConstOp(2.0) * sqrt(expr.arg)) {}

  template <typename Data>
  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return derivative(x1, x2, x3, x4);
  }
};

// Derivative of a constant
template <int Argument>
struct Derivative<Argument, ConstOp >
{
  typedef ConstOp arg_type;
  typedef ZeroOp result_type;
  result_type derivative;

  HOST_DEVICE Derivative() : derivative() {}

  template <typename Data>
  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return derivative(x1, x2, x3, x4);
  }
};

template <int Argument>
struct Derivative<Argument, double>
{
  typedef double arg_type;
  typedef ZeroOp result_type;
  result_type derivative;

  HOST_DEVICE Derivative() : derivative() {}

  template <typename Data>
  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return derivative(x1, x2, x3, x4);
  }
};

template <int Argument>
struct Derivative<Argument, ZeroOp>
{
  typedef ZeroOp arg_type;
  typedef ZeroOp result_type;
  result_type derivative;

  HOST_DEVICE Derivative() : derivative() {}

  template <typename Data>
  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return derivative(x1, x2, x3, x4);
  }
};

// Derivative of an independent variable, producing Kronicker delta
template <int Argument, int var, typename Data>
struct Derivative<Argument, Var<var, Data> >
{
  typedef Var<var, Data> arg_type;
  typedef ZeroOp result_type;
  result_type derivative;

  HOST_DEVICE Derivative() {
    static_assert(Argument != var, "Template matching error, Var and argument should be same");
    derivative = ZeroOp();
  }

  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return derivative(x1, x2, x3, x4);
  }
};

template <int var, typename Data>
struct Derivative<var, Var<var, Data> >
{
  typedef Var<var, Data> arg_type;
  typedef ConstOp result_type;
  result_type derivative;

  HOST_DEVICE Derivative() {
    derivative = ConstOp(1.0);
  }

  HD_INLINE Data operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0) {
    return derivative(x1, x2, x3, x4);
  }
};

// Multiple derivatives
// template <int Argument, int times, typename Expr>
// struct Derivative_n
// {
//     typedef Expr arg_type;
//     typedef typename Derivative<Argument,
//                                 typename Derivative_n<Argument, times - 1, Expr>::result_type>::result_type result_type;
//     result_type derivative;

//     HOST_DEVICE Derivative_n(arg_type expr) : derivative(D<Argument>(D_n<Argument, times - 1>(expr))) {
//     }

//     HD_INLINE double operator() (double x1, double x2 = 0.0, double x3 = 0.0, double x4 = 0.0) {
//         return derivative(x1, x2, x3, x4);
//     }
// };

// template <int Argument, typename Expr>
// struct Derivative_n<Argument, 1, Expr>
// {
//     typedef Expr arg_type;
//     typedef typename Derivative<Argument, Expr>::result_type result_type;
//     result_type derivative;

//     HOST_DEVICE Derivative_n(arg_type expr) : derivative(D<Argument>(expr)) {}

//     HD_INLINE double operator() (double x1, double x2 = 0.0, double x3 = 0.0, double x4 = 0.0) {
//         return derivative(x1, x2, x3, x4);
//     }
// };

// template <int Argument, typename Expr>
// struct Derivative_n<Argument, 0, Expr>
// {
//     typedef Expr arg_type;
//     typedef Expr result_type;
//     result_type derivative;

//     HOST_DEVICE Derivative_n(arg_type expr) : derivative(expr) {}

//     HD_INLINE double operator() (double x1, double x2 = 0.0, double x3 = 0.0, double x4 = 0.0) {
//         return derivative(x1, x2, x3, x4);
//     }
// };


// template <int Arg1, int n1, int Arg2, int n2, typename Expr>
// struct Derivative_2n
// {
//     typedef Expr arg_type;
//     typedef typename Derivative<Arg1,
//                                 typename Derivative_2n<Arg1, n1 - 1, Arg2, n2, Expr>::result_type>::result_type result_type;
//     result_type derivative;

//     HOST_DEVICE Derivative_2n(arg_type expr) : derivative(D<Arg1>(D_nm<Arg1, n1 - 1, Arg2, n2>(expr))) {}

//     HD_INLINE double operator() (double x1, double x2 = 0.0, double x3 = 0.0, double x4 = 0.0) {
//         return derivative(x1, x2, x3, x4);
//     }
// };

// template <int Arg1, int Arg2, int n2, typename Expr>
// struct Derivative_2n<Arg1, 1, Arg2, n2, Expr>
// {
//     typedef Expr arg_type;
//     typedef typename Derivative<Arg1,
//                                 typename Derivative_n<Arg2, n2, Expr>::result_type>::result_type result_type;
//     result_type derivative;

//     HOST_DEVICE Derivative_2n(arg_type expr) : derivative(D<Arg1>(D_n<Arg2, n2>(expr))) {}

//     HD_INLINE double operator() (double x1, double x2 = 0.0, double x3 = 0.0, double x4 = 0.0) {
//         return derivative(x1, x2, x3, x4);
//     }
// };

// template <int Arg1, int Arg2, int n2, typename Expr>
// struct Derivative_2n<Arg1, 0, Arg2, n2, Expr>
// {
//     typedef Expr arg_type;
//     typedef typename Derivative_n<Arg2, n2, Expr>::result_type result_type;
//     result_type derivative;

//     HOST_DEVICE Derivative_2n(arg_type expr) : derivative(D_n<Arg2, n2>(expr)) {}

//     HD_INLINE double operator() (double x1, double x2 = 0.0, double x3 = 0.0, double x4 = 0.0) {
//         return derivative(x1, x2, x3, x4);
//     }
// };

////////////////////////////////////////////////////////////////////////////////
///  Functions
////////////////////////////////////////////////////////////////////////////////
template <int Argument, typename Op, typename Left, typename Right>
HD_INLINE
typename Derivative<Argument, BinaryOp<Op, Left, Right> >::result_type
D(const BinaryOp<Op, Left, Right>& expr) {
  // printf("Constructed derivative: %s\n", __PRETTY_FUNCTION__);
  return Derivative<Argument, BinaryOp<Op, Left, Right> >(expr.left, expr.right).derivative;
}

template <int Argument,typename Op, typename Right>
HD_INLINE
typename Derivative<Argument, BinaryOp<Op, double, Right> >::result_type
D(const BinaryOp<Op, double, Right>& expr) {
  return Derivative<Argument, BinaryOp<Op, double, Right> >(expr.left.val, expr.right).derivative;
}

template <int Argument,typename Op, typename Left>
HD_INLINE
typename Derivative<Argument, BinaryOp<Op, Left, double> >::result_type
D(const BinaryOp<Op, Left, double>& expr) {
  return Derivative<Argument, BinaryOp<Op, Left, double> >(expr.left, expr.right.val).derivative;
}

template <int Argument,typename Op, typename Expr>
HD_INLINE
typename Derivative<Argument, UnaryOp<Op, Expr> >::result_type
D(const UnaryOp<Op, Expr>& expr) {
  return Derivative<Argument, UnaryOp<Op, Expr> >(expr.arg).derivative;
}

// template <int Argument, int n>
// HD_INLINE
// typename Derivative<Argument, Var<n, double> >::result_type
// D(const Var<n, double>& var) {
//   return Derivative<Argument, Var<n, double> >().derivative;
// }

template <int Argument, int n, typename Data>
HD_INLINE
typename Derivative<Argument, Var<n, Data> >::result_type
D(const Var<n, Data>& var) {
  return Derivative<Argument, Var<n, Data> >().derivative;
}

template <int Argument>
HD_INLINE
typename Derivative<Argument, double>::result_type
D(const double& val) {
  return Derivative<Argument, double>().derivative;
}

template <int Argument>
HD_INLINE
typename Derivative<Argument, ConstOp>::result_type
D(const ConstOp& val) {
  return Derivative<Argument, ConstOp>().derivative;
}

template <int Argument>
HD_INLINE
typename Derivative<Argument, ZeroOp>::result_type
D(const ZeroOp& val) {
  return Derivative<Argument, ZeroOp>().derivative;
}

template <int Argument, int Arg2, typename Expr>
HD_INLINE
typename Derivative<Argument, typename Derivative<Arg2, Expr>::result_type>::result_type
D(const Derivative<Arg2, Expr>& deriv) {
  return D<Argument>(deriv.derivative);
}

// template <int Argument, int times, typename Expr>
// HD_INLINE
// typename Derivative_n<Argument, times, Expr>::result_type
// D_n_imp(const Expr& expr, int2type<times> I) {
//     return D<Argument>(D_n_imp<Argument>(expr, int2type<times - 1>()));
// }

// template <int Argument, typename Expr>
// HD_INLINE
// typename Derivative<Argument, Expr>::result_type
// D_n_imp(const Expr& expr, int2type<1> I) {
//     return D<Argument>(expr);
// }

// template <int Argument, typename Expr>
// HD_INLINE
// Expr
// D_n_imp(const Expr& expr, int2type<0> I) {
//     return expr;
// }

// template <int Argument, int times, typename Expr>
// HD_INLINE
// typename Derivative_n<Argument, times, Expr>::result_type
// D_n(const Expr& expr) {
//     return D_n_imp<Argument>(expr, int2type<times>());
// }

// template <int Arg1, int n1, int Arg2, int n2, typename Expr>
// HD_INLINE
// typename Derivative_2n<Arg1, n1, Arg2, n2, Expr>::result_type
// D_nm_imp(const Expr& expr, int2type<n1> I, int2type<n2> J) {
//     return D<Arg1>(D_nm_imp<Arg1, n1 - 1, Arg2, n2>(expr, int2type<n1 - 1>(), int2type<n2>()));
// }

// template <int Arg1, int Arg2, int n2, typename Expr>
// HD_INLINE
// typename Derivative<Arg1, typename Derivative_n<Arg2, n2, Expr>::result_type>::result_type
// D_nm_imp(const Expr& expr, int2type<1> I, int2type<n2> J) {
//     return D<Arg1>(D_n<Arg2, n2>(expr));
// }

// template <int Arg1, int n1, int Arg2, int n2, typename Expr>
// HD_INLINE
// typename Derivative_2n<Arg1, n1, Arg2, n2, Expr>::result_type
// D_nm(const Expr& expr) {
//     return D_nm_imp<Arg1, n1, Arg2, n2>(expr, int2type<n1>(), int2type<n2>());
// }


}

#endif   // ----- #ifndef _DERIVATIVE_H_  -----
