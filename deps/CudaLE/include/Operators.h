/////////////////////////////////////////////////////////////////////////////////////////
///
///           \file  Operators.h
///
/// __Description__:     Declares the base operators for the expression template
///
/// __Last modified__:   <2014-07-09 01:46:36 alex>\n
/// __Version__:         1.0\n
/// __Author__:          Alex Chen, fizban007@gmail.com\n
/// __Organization__:    Columbia University
///
/////////////////////////////////////////////////////////////////////////////////////////

#ifndef  _OPERATORS_H_
#define  _OPERATORS_H_

#include <type_traits>
#include "cudaControl.h"
#include "helper.h"

namespace CudaLE {

struct ZeroOp
{
  typedef ZeroOp type;
  HOST_DEVICE ZeroOp() {}
  template <typename Any>
  HOST_DEVICE ZeroOp(const Any& any) {}
  template <typename Any1, typename Any2>
  HOST_DEVICE ZeroOp(const Any1& any1, const Any2& any2) {}

  template <typename Data>
  HD_INLINE double operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0, const Data& x5 = 0.0, const Data& x6 = 0.0) const {
    return 0.0;
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value, bool>::type
  operator== (const T& obj) const {
    return true;
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value == false, bool>::type
  operator== (const T& obj) const {
    return false;
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value, bool>::type
  operator!= (const T& obj) const {
    return false;
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value == false, bool>::type
  operator!= (const T& obj) const {
    return true;
  }

  HD_INLINE void print() const {
    helper::print("Zero");
  }
};

struct OneOp
{
  typedef OneOp type;

  HOST_DEVICE OneOp() {}
  template <typename Any>
  HOST_DEVICE OneOp(const Any& any) {}
  template <typename Any1, typename Any2>
  HOST_DEVICE OneOp(const Any1& any1, const Any2& any2) {}

  template <typename Data>
  HD_INLINE double operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0, const Data& x5 = 0.0, const Data& x6 = 0.0) const {
    return 1.0;
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value, bool>::type
  operator== (const T& obj) const {
    return true;
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value == false, bool>::type
  operator== (const T& obj) const {
    return false;
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value, bool>::type
  operator!= (const T& obj) const {
    return false;
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value == false, bool>::type
  operator!= (const T& obj) const {
    return true;
  }

  HD_INLINE void print() const {
    helper::print(1.0);
  }
};

struct ConstOp
{
  typedef ConstOp type;

  double val;
  HOST_DEVICE ConstOp(double v) : val(v) {}
  HOST_DEVICE ConstOp(const ConstOp& op) : val(op.val) {}
  HOST_DEVICE ConstOp() : val(0.0) {}

  // Allow ConstOp to be converted to double automatically
  HD_INLINE operator double() const { return val; }

  template <typename Data>
  HD_INLINE double operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0, const Data& x5 = 0.0, const Data& x6 = 0.0) const {
    return val;
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value, bool>::type
  operator== (const T& obj) const {
    return obj.val == val;
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value == false, bool>::type
  operator== (const T& obj) const {
    return false;
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value, bool>::type
  operator!= (const T& obj) const {
    return (val != obj.val);
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value == false, bool>::type
  operator!= (const T& obj) const {
    return true;
  }

  HD_INLINE void print() const {
    helper::print(val);
  }
};

template <typename Op, typename Left, typename Right>
struct BinaryOp
{
  Left left;
  Right right;
  typedef BinaryOp<Op, Left, Right> type;

  HOST_DEVICE BinaryOp(Left t1, Right t2) : left(t1), right(t2) {}
  HOST_DEVICE BinaryOp(const type& op) : left(op.left), right(op.right) {}
  HOST_DEVICE BinaryOp() {}

  template <typename Data>
  HD_INLINE auto operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0, const Data& x5 = 0.0, const Data& x6 = 0.0) const {
    return Op::apply(left(x1, x2, x3, x4, x5, x6), right(x1, x2, x3, x4, x5, x6));
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value, bool>::type
  operator== (const T& obj) const {
    return ((left == obj.left) && (right == obj.right));
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value == false, bool>::type
  operator== (const T& obj) const {
    return false;
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value, bool>::type
  operator!= (const T& obj) const {
    return (!operator==(obj));
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value == false, bool>::type
  operator!= (const T& obj) const {
    return true;
  }

  HD_INLINE void print() const {
    helper::print("(");
    left.print();
    Op::print();
    right.print();
    helper::print(")");
  }
};

template <typename Op, typename Left>
struct BinaryOp<Op, Left, double>
{
  Left left;
  ConstOp right;
  typedef BinaryOp<Op, Left, double> type;

  HOST_DEVICE BinaryOp(Left t1, double t2) : left(t1), right(ConstOp(t2)) {}
  HOST_DEVICE BinaryOp(Left t1, ConstOp t2) : left(t1), right(t2) {}
  HOST_DEVICE BinaryOp(const type& op) : left(op.left), right(op.right) {}
  HOST_DEVICE BinaryOp() {}

  template<typename Data>
  HD_INLINE auto operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0, const Data& x5 = 0.0, const Data& x6 = 0.0) const {
    return Op::apply(left(x1, x2, x3, x4, x5, x6), right(x1, x2, x3, x4, x5, x6));
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value, bool>::type
  operator== (const T& obj) const {
    return ((left == obj.left) && (right == obj.right));
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value == false, bool>::type
  operator== (const T& obj) const {
    return false;
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value, bool>::type
  operator!= (const T& obj) const {
    return (!operator==(obj));
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value == false, bool>::type
  operator!= (const T& obj) const {
    return true;
  }

  HD_INLINE void print() const {
    helper::print("(");
    left.print();
    Op::print();
    right.print();
    helper::print(")");
  }
};

template <typename Op, typename Right>
struct BinaryOp<Op, double, Right>
{
  ConstOp left;
  Right right;
  typedef BinaryOp<Op, double, Right> type;

  HOST_DEVICE BinaryOp(double t1, Right t2) : left(ConstOp(t1)), right(t2) {}
  HOST_DEVICE BinaryOp(ConstOp t1, Right t2) : left(t1), right(t2) {}
  HOST_DEVICE BinaryOp(const type& op) : left(op.left), right(op.right) {}
  HOST_DEVICE BinaryOp() {}

  template <typename Data>
  HD_INLINE auto operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0, const Data& x5 = 0.0, const Data& x6 = 0.0) const {
    return Op::apply(left(x1, x2, x3, x4, x5, x6), right(x1, x2, x3, x4, x5, x6));
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value, bool>::type
  operator== (const T& obj) const {
    return ((left == obj.left) && (right == obj.right));
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value == false, bool>::type
  operator== (const T& obj) const {
    return false;
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value, bool>::type
  operator!= (const T& obj) const {
    return (!operator==(obj));
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value == false, bool>::type
  operator!= (const T& obj) const {
    return true;
  }

  HD_INLINE void print() const {
    helper::print("(");
    left.print();
    Op::print();
    right.print();
    helper::print(")");
  }
};

template <typename Op>
struct BinaryOp<Op, double, double>
{
  ConstOp left;
  ConstOp right;
  typedef BinaryOp<Op, double, double> type;

  HOST_DEVICE BinaryOp(double t1, double t2) : left(ConstOp(t1)), right(ConstOp(t2)) {}
  HOST_DEVICE BinaryOp(ConstOp t1, ConstOp t2) : left(t1), right(t2) {}
  HOST_DEVICE BinaryOp(const type& op) : left(op.left), right(op.right) {}
  HOST_DEVICE BinaryOp() {}

  template<typename Data>
  HD_INLINE auto operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0, const Data& x5 = 0.0, const Data& x6 = 0.0) const {
    return Op::apply(left(x1, x2, x3, x4, x5, x6), right(x1, x2, x3, x4, x5, x6));
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value, bool>::type
  operator== (const T& obj) const {
    return ((left == obj.left) && (right == obj.right));
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value == false, bool>::type
  operator== (const T& obj) const {
    return false;
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value, bool>::type
  operator!= (const T& obj) const {
    return (!operator==(obj));
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value == false, bool>::type
  operator!= (const T& obj) const {
    return true;
  }

  HD_INLINE void print() const {
    helper::print("(");
    left.print();
    Op::print();
    right.print();
    helper::print(")");
  }
};

template <typename Op, typename Arg>
struct UnaryOp
{
  Arg arg;
  typedef UnaryOp<Op, Arg> type;

  HOST_DEVICE UnaryOp(Arg t1) : arg(t1) {}
  HOST_DEVICE UnaryOp(const type& op) : arg(op.arg) {}
  HOST_DEVICE UnaryOp() {}

  template <typename Data>
  HD_INLINE auto operator() (const Data& x1, const Data& x2 = 0.0, const Data& x3 = 0.0, const Data& x4 = 0.0, const Data& x5 = 0.0, const Data& x6 = 0.0) const {
    return Op::apply(arg(x1, x2, x3, x4, x5, x6));
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value, bool>::type
  operator== (const T& obj) const {
    return (arg == obj.arg);
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value == false, bool>::type
  operator== (const T& obj) const {
    return false;
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value, bool>::type
  operator!= (const T& obj) const {
    return (!operator==(obj));
  }

  template <typename T>
  HD_INLINE
  typename std::enable_if<std::is_same<T, type>::value == false, bool>::type
  operator!= (const T& obj) const {
    return true;
  }

  HD_INLINE void print() const {
    Op::print();
    helper::print("(");
    arg.print();
    helper::print(")");
  }
};

// template <typename Op, typename Arg>
// struct UnaryOp<Op, Arg*>
// {
//     Arg* arg;
//     typedef UnaryOp<Op, Arg*> type;

//     HOST_DEVICE UnaryOp(Arg* t1) : arg(t1) {}
//     HOST_DEVICE UnaryOp() {}

//     HD_INLINE double operator() (double x1, double x2 = 0.0, double x3 = 0.0, double x4 = 0.0, const Data& x5 = 0.0, const Data& x6 = 0.0) {
//         return Op::apply((*arg)(x1, x2, x3, x4, x5, x6));
//     }
// };
}

#endif   // ----- #ifndef _OPERATORS_H_  -----
