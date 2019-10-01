//==================================================================================
// Module       : gtstat.hpp
// Date         : 10/1/19 (DLR)
// Description  : Encapsulates the methods and data associated with
//                GeoFLOW statistics methods.
// Copyright    : Copyright 2019, Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================

#if !defined(_GTSTAT_DECL_HPP)
#define _GTSTAT_DECL_HPP

#include <cstdlib>
#include <limits>
#include <iostream>
#include <vector>
#include "gtypes.h"
#include "gtvector.hpp"
#include "gtmatrix.hpp"
#include "gcomm.hpp"



template <class T> class GTStat
{
  assert(std::is_floating_point<T>::value && "Requires floating point template parameter");

  public:

           GTStat<T>() = delete;
           GTStat<T>(GSIZET nbins, GC_COMM icomm=GC_COMM_WORLD);
          ~GTStat<T>() = default;
    
    void   set_nbins(GSIZE nbins)
             {nbins_ = nbin;)             // Set no. bins
    void   get_nbins(GSIZE nbin)
             {return nbins_ );            // Get no. bins
    void   dopdf1d(GTVector<T> u, GBOOL ifixdr, T &fmin, T &fmax, GDOOL dolog, GTVector<T> &pdf);
    void   dopdf1d(GTVector<T> u, GBOOL ifixdr, T &fmin, T &fmax, GDOOL dolog, const GString &fname); 
    

  private:

    GINT             myrank_;   // task's rank
    GSIZET           nbins_;    // no. bins
    GSIZET           nkeep_;    // no. stoch var indices kept
    GTVector<GSIZET> ikeep_;    // indirection indices to valid values
    GTVector     <T> lpdf_ ;    // local pdf
    GTVector     <T> gpdf_ ;    // globally-reduced pdf

};


#include "gtstat.ipp"

#endif

