#include<iostream>
#include<math.h>
#include "CudaLE.h"

using namespace CudaLE;

// const CoordinateType type = COORD_CYLINDRICAL;

// typedef UnaryOp<Id, double> Var;
template <typename T>
struct ScaleFunc
{
    T h;
    ScaleFunc (T func) : h(func) {}
    double operator() (double x1, double x2) {
        return h(x1, x2);
    }
}; // ----- end of class Scales -----

template <typename T>
ScaleFunc<T> MakeScaleFunc(T func) {
    return ScaleFunc<T>(func);
}

template <typename T1, typename T2, typename T3>
struct Scales
{
    T1 h1;
    T2 h2;
    T3 h3;
    Scales(T1 f1, T2 f2, T3 f3) : h1(f1), h2(f2), h3(f3) {}
    template <int i>
    double h (double x1, double x2) {
        if (1 == i) return h1(x1, x2);
        else if (2 == i) return h2(x1, x2);
        else if (3 == i) return h3(x1, x2);
        else return 1.0;
    }
}; // ----- end of class Scales -----


template <typename T1, typename T2, typename T3>
Scales<T1, T2, T3> MakeScales(T1 f1, T2 f2, T3 f3) {
    return Scales<T1, T2, T3>(f1, f2, f3);
}


static auto myScalesCart = MakeScales(ConstOp(1.0), ConstOp(1.0), ConstOp(1.0));

int main () 
{
    // Var _x(1), _y(2);
    // (_sin(_1) + _2)::op_type myF = MakeFunc(_sin(_1) + _2);
    auto myF = _sin(_1) + _2;
    auto myG = _1 * 4.0 - 2.0 *_2;
    auto myH = _1 / 4.0 - 2.0 /_2;
    auto myComb = myF * myG;
    auto myCart = ConstOp(1.0);
    auto myScalesSpherical = MakeScales(ConstOp(1.0), _1, _1 * _sin(_2));
    

    std::cout << myScalesCart.h<2>(2.0,324.0) << std::endl;
    std::cout << myScalesSpherical.h<3>(23.0,2.0) << std::endl;
    std::cout << (myScalesSpherical.h3 * myScalesCart.h2)(23.1, 3.0) << std::endl;

    auto myScale1 = MakeScaleFunc(myF);
    // std::cout << sin(2.0) + 3.0 << std::endl;
    // std::cout << (_log(_1) * _exp(_2) + _3)(2.0, -1, 1.0) << std::endl;
    std::cout << myF(2.0, 2.0) << std::endl;
    std::cout << myG(2.0, 2.0) << std::endl;
    std::cout << myComb(2.0, 2.0) << std::endl;
    std::cout << myScale1(2.0, 2.0) << std::endl;
    std::cout << myCart(1.0, 2.0, 10.0) << std::endl;
    return 0;
}
