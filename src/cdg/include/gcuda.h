//==================================================================================
// Module       : gcuda.hpp
// Date         : 12/3/20 (DLR)
// Description  : Basic types, defs for CUDA interfaces
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GCUDA_HPP)
#define _GCUDA_HPP

#include "gtvector.hpp"

#if defined(_G_USE_CUDA)

  #define GCuStream    cudaStream_t
  #define GBlasHandle  cublasHandle_t
  #define GBlasStatus  cublasStatus_t
  #define GBlasOp      cublasOperation_t
  #define GBlasError   cudaError_t

#else

  #define GCuStream    void*
  #define GBlasHandle  void*
  #define GBlasStatus  void*
  #define GBlasOp      int
  #define GBlasError   int

#endif

struct cuMatBlockDat {
  GTVector        <GINT> ibblk;  // for each stream, starting blk index
  GTVector        <GINT> ieblk;  // for each stream, ending blk index
  GTVector<GCuStream>  pStream;  // stream pointers
};


#endif // !defined(_GTYPES_HPP)

