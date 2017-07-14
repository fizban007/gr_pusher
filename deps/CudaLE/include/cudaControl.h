/////////////////////////////////////////////////////////////////////////////////////////
///
///           \file  cudaControl.h
///
/// __Description__:     Defines some useful macros
///
/// __Version__:         1.0\n
/// __Author__:          Alex Chen, fizban007@gmail.com\n
/// __Organization__:    Columbia University
///
/////////////////////////////////////////////////////////////////////////////////////////

#ifndef  _CUDACONTROL_H_
#define  _CUDACONTROL_H_

#ifdef __CUDACC__
  #define HOST_DEVICE __host__ __device__
  #define HD_INLINE __host__ __device__ __forceinline__
  #define WITH_CUDA_ENABLED
#else
  #define HOST_DEVICE
  #define HD_INLINE inline
  #ifdef WITH_CUDA_ENABLED
    #undef WITH_CUDA_ENABLED
  #endif
#endif   // ----- #ifdef __CUDACC__ -----

#endif   // ----- #ifndef _CUDACONTROL_H_  -----
