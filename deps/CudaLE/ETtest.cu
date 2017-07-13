
#include<math.h>
    
#ifdef __CUDACC__
#define HOST_DEVICE __host__ __device__
#else
#define HOST_DEVICE
#endif
    
struct Id {};
    
template <typename Op, typename Left, typename Right>
struct BinaryOp
{
    Left left;
    Right right;
    HOST_DEVICE BinaryOp(Left t1, Right t2) : left(t1), right(t2) {}
    
    HOST_DEVICE double operator() (double x1, double x2 = 0.0, double x3 = 0.0, double x4 = 0.0) {
        return Op::apply(left(x1, x2, x3), right(x1, x2, x3));
    }
};
    
template <typename Op, typename Arg>
struct UnaryOp
{
    Arg arg;
    HOST_DEVICE UnaryOp(Arg t1) : arg(t1) {}
    
    HOST_DEVICE double operator() (double x1, double x2 = 0.0, double x3 = 0.0, double x4 = 0.0) {
        return Op::apply(arg(x1, x2, x3));
    }
};
    
template <int argnum>
struct Var
{
    HOST_DEVICE Var() {}
    HOST_DEVICE double operator() (double x1, double x2 = 0.0, double x3 = 0.0, double x4 = 0.0) {
        if (1 == argnum) return x1;
        else if (2 == argnum) return x2;
        else return x3;
    }
};
    
struct Sin
{
    HOST_DEVICE static double apply(double x) {
        return sin(x);
    }
};
    
struct Plus
{
    HOST_DEVICE static double apply(double a, double b) {
        return a + b;
    }
};
    
template <typename Left, typename Right>
BinaryOp<Plus, Left, Right> operator+ (Left lhs, Right rhs) {
    return BinaryOp<Plus, Left, Right>(lhs, rhs);
}
    
template <typename Arg>
UnaryOp<Sin, Arg> _sin(Arg arg) {
    return UnaryOp<Sin, Arg>(arg);
}
    
template <class T>
__global__ void test(T func, double x, double y, double z = 0.0) {
    printf("%e\n", func(x, y));
}

Var<1> _x;
Var<2> _y;
    
int main () 
{
    test<<<1, 1>>>(_sin(_x) + _y, 1.0, 2.0);
    cudaDeviceSynchronize();  // Needed or the host will return before kernel is finished
    return 0;
}
