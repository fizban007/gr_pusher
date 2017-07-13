#include<iostream>
#include<math.h>
#include<typeinfo>
#include<boost/detail/sp_typeinfo.hpp>
#include "CudaLE.h"

using namespace CudaLE;

#define CURRENT_FUNCTION __PRETTY_FUNCTION__

// const Coordinate Type type = COORD_CYLINDRICAL;
// struct ScaleFuncBase
// {
//     HOST_DEVICE virtual double operator() (double x1, double x2 = 0.0, double x3 = 0.0, double x4 = 0.0) = 0;
//     HOST_DEVICE virtual ~ScaleFuncBase() {}
// };

// // typedef UnaryOp<Id, double> Var;
// template <typename Func>
// struct ScaleFunc : public ScaleFuncBase
// {
//     Func h;
//     HOST_DEVICE ScaleFunc (Func func) : h(func) {}
//     HOST_DEVICE virtual double operator() (double x1, double x2 = 0.0, double x3 = 0.0, double x4 = 0.0) {
//         return h(x1, x2);
//     }
//     HOST_DEVICE virtual ~ScaleFunc() {}
// }; // ----- end of class Scales -----

// template <typename F>
// struct Finvoke
// {
//     HOST_DEVICE static double invoke (F func, double x1, double x2, double x3, double x4 = 0.0) {
//         F f(func);
//         return f(x1, x2, x3);
//     }
// };

struct FunctionBase
{
    HOST_DEVICE virtual double operator() (double x1, double x2 = 0.0, double x3 = 0.0, double x4 = 0.0) = 0;
    HOST_DEVICE virtual ~FunctionBase() {}
    // template<typename F>
    // F getFunc(FunctionBase*);
};

template <typename F>
struct Function : public FunctionBase
{
    F f;
    HOST_DEVICE Function(F func) : f(func) {}
    HOST_DEVICE virtual double operator() (double x1, double x2 = 0.0, double x3 = 0.0, double x4 = 0.0) {
        return f(x1, x2, x3);
    }
    HOST_DEVICE virtual ~Function() {}
};

// template <typename F>
// struct func_wrapper : public F
// {
//     func_wrapper(F f) : F(f) {}
// };

// struct function_buffer
// {
//     void * f;

//     struct type_t {
//         // const boost::detail::sp_typeinfo* type;
//         std::type_info* type;
//     } type;
// };

// template<typename Functor>
// struct functor_manager
// {
//     HOST_DEVICE static inline void
//     manage(const function_buffer& in_buffer) {
//         Functor* f = static_cast<Functor*>(in_buffer.f);
//         delete f;
//         // in_buffer.f = 0;
//     }
// };

// template <typename F>
// struct Invoker
// {
//     HOST_DEVICE static double invoke(function_buffer& functor, double x1, double x2, double x3, double x4 = 0.0) {
//         F* f;
//         f = reinterpret_cast<F*>(functor.f);
//         return (*f)(x1, x2, x3);
//     }
// };

// struct Get_Invoker
// {
//     template<typename F>
//     struct apply
//     {
//         typedef Invoker<F> invoker_type;
//         typedef functor_manager<F> manage_type;
//     };
// };

// struct vtable_base
// {
//     void (*manager)(const function_buffer& in_buffer);
// };

// struct vtable
// {
//     typedef double (*invoker_type)(function_buffer&,
//                                    double, double, double);

//     void (*manager)(const function_buffer& in_buffer);
    
//     template<typename F>
//     HOST_DEVICE void
//     assign_functor(F f, function_buffer& functor) const {
//         functor.f = new F(f);
//     }

//     void
//     clean (function_buffer& functor) const {
//     }

//     invoker_type invoker;
// };

// struct FunctionTest
// {
//     function_buffer functor;
//     vtable* vtable_;

//     template <typename F>
//     HOST_DEVICE FunctionTest(F func) {
//         assign(func);
//     }

//     HOST_DEVICE double operator() (double x1, double x2 = 0.0, double x3 = 0.0, double x4 = 0.0) {
//         return vtable_ -> invoker(functor, x1, x2, x3);
//     }

//     template <typename F>
//     HOST_DEVICE void assign(F f) {
//         typedef typename Get_Invoker::template apply<F> handler_type;
//         typedef typename handler_type::invoker_type invoker_type;
//         typedef typename handler_type::manage_type manager_type;

//         // static const vtable stored_vtable = { &invoker_type::invoke, &manager_type::manage };
//         vtable stored_vtable;
//         stored_vtable.invoker = &invoker_type::invoke;
//         stored_vtable.manager = &manager_type::manage;
//         stored_vtable.assign_functor(f, functor);
//         // functor.type.type = &typeid(F);
//         size_t value = reinterpret_cast<size_t>(&stored_vtable);
//         vtable_ = reinterpret_cast<vtable*>(value);
//     }
// };

// template <typename T>
// ScaleFuncBase* MakeScaleFunc(T func) {
//     return (new ScaleFunc<T>(func));
// }

template <typename T>
FunctionBase* MakeFunction(T func) {
    FunctionBase* newFunc = new Function<T>(func);
    // newFunc -> invoker_type = &Finvoke<T>::invoke;
    return newFunc;
}

// template<class T> struct sp_typeid_ {
//     static char const * name() {
//         return CURRENT_FUNCTION;
//     }
// };

// typedef UnaryOp<Id, double> Var;
template <typename F>
__global__ void test (F func, double x) {
    // printf("%s\n", typeid(int).name());
    // printf("%s\n", sp_typeid_::name());
    printf("%e\n", func(x));
}

int main () 
{
    // FunctionBase* myF = MakeFunction(_sin(_1) + _2);
    FunctionBase& myG = *MakeFunction(_sin(_2) * _exp(_1) + _2 + _cos(_3));
    // typeof(*myG) mytest;
    // Function func(_sin(_1));
    // std::cout << "Testing vtable" << func(2.0) << std::endl;
    // std::cout << typeid(_sin(_1)).name() << std::endl;
    // Var _x(1), _y(2);
    // (_sin(_1) + _2)::op_type myF = MakeFunc(_sin(_1) + _2);
    // auto myF = _sin(_1) + _2;
    // std::cout << sin(2.0) + 3.0 << std::endl;
    // std::cout << (_log(_1) * _exp(_2) + _3)(2.0, -1, 1.0) << std::endl;
    // std::cout << _sin(_pow<2>(_1) + _2)(2.0, 3.0) << std::endl;
    // std::cout << (myG -> getFunc(myG))(2.0, 3.0) << std::endl;
    std::cout << myG(1.0, 2.0, 3.0) << std::endl;
    delete &myG;
    // std::cout << (myF * _sin(_2))(1.0, 2.0) << std::endl;
    // std::cout << (_exp(myG))(1.0, 2.0) << std::endl;
    // std::cout << sizeof(myG) << std::endl;
    // test<<<1, 1>>>(_sin(_1), 2.0);
    // test<<<1, 1>>>((* dynamic_cast<Function<typeof(_sin(_2) + _1)>*>(myG)).f, 2.0);
    cudaDeviceSynchronize();
    // delete myF;
    // delete myG;
    return 0;
}
