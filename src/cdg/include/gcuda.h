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

  #define GCUSTREAM cudaStream_t

#else

  #define GCUSTREAM void*

#endif

struct cuMatBlockDat {
  GTVector        <GINT> ibblk;  // for each stream, starting blk index
  GTVector        <GINT> ieblk;  // for each stream, ending blk index
  GTVector<GCUSTREAM>  pStream;  // stream pointers
};


#endif // !defined(_GTYPES_HPP)

