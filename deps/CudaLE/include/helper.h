/////////////////////////////////////////////////////////////////////////////////////////
///
///           \file  helper.h
///
/// __Description__:     Defines helper functions like print().
///
/// __Last modified__:   <>\n
/// __Version__:         1.0\n
/// __Author__:          Alex Chen, fizban007@gmail.com\n
/// __Organization__:    Columbia University
///
/////////////////////////////////////////////////////////////////////////////////////////


#ifndef  _HELPER_H_
#define  _HELPER_H_

#include <iostream>
#include "cudaControl.h"

namespace CudaLE {

namespace helper {

HD_INLINE void print(const char* str) {
#ifdef WITH_CUDA_ENABLED
  printf("%s", str);
#else
  std::cout << str;
#endif
}

HD_INLINE void print(double d) {
#ifdef WITH_CUDA_ENABLED
  printf("%f", d);
#else
  std::cout << d;
#endif
}

template <typename T>
HD_INLINE void println(const T& t) {
  t.print();
  print("\n");
}

}

}

#endif   // _HELPER_H_
