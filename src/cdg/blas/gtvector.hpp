//==================================================================================
// Module       : gtvector_decl.hpp
// Date         : 1/1/18 (DLR)
// Description  : Encapsulates the methods and data associated with
//                a template vector object, with contiguous memory
//                for basic types.
//                NOTE: Partial support for generalized 'index' that
//                      yields, start, stop, stride, pad indices.
//                The implementation file, gtvector.ipp, is not included
//                in this file, just the declarations.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

#if !defined(_GTVECTOR_DECL_HPP)
#define _GTVECTOR_DECL_HPP

#include "gtypes.h"
#include <cstdlib>
#include <limits>
#include <iostream>
#include "gindex.hpp"
#include "cff_blas.h"

#if !defined(_G_VEC_CACHE_SIZE)
  # define _G_VEC_CACHE_SIZE 16
#endif



template <class T> class GTVector
{
  public:

    GTVector<T>();
    GTVector<T>(GSIZET n);
    GTVector<T>(GIndex &gin);
    GTVector<T>(GTVector<T> &obj);
    GTVector<T>(const GTVector<T> &obj);
    GTVector<T>(T *, GSIZET n, GSIZET istride=1);
    GTVector<T>(T *, GSIZET n, GSIZET istride, GBOOL bmanaged=TRUE);
   ~GTVector<T>();
    
    T *data();
    GSIZET size() const;        // Return used buffer size
    GSIZET capacity() const;    // Return total available buffer size
    void   reserve(GSIZET n);   // Set capacity before having to initialize
    void   resize(GSIZET n);    // Resize data buffer
    void   resize(GIndex &);    // Resize data buffer with gindex
    void   resizem(GSIZET n);   // Resize data buffer only if n > current
    void   clear();             // Set capacity to 0
    void   push_back(T const &);// Push new element to end of data buffer
    void   bconstdata(GBOOL);   // If data const, all access access 1 element only

    void range(GSIZET ibeg, GSIZET end);    // Set range of vector within capacity
    GIndex &getIndex() ;        // Return generalized index member

    #pragma acc routine vector
    GTVector<T>       &operator=(const GTVector<T> &b);
    #pragma acc routine vector
    void               operator=(T b);
    #pragma acc routine vector
    void               pointProd(const GTVector<T> &fact, GTVector<T> &ret);
    #pragma acc routine vector
    void               pointProd(const GTVector<T> &);
    #pragma acc routine vector
    void               constProd(const T, GTVector<T> &ret);
    #pragma acc routine vector
    void transpose(GSIZET n);

  
    #pragma acc routine vector
    void               set(T b);
    #pragma acc routine vector
    void               set(T *b, GSIZET n);
    // Device//accelerator data methods:
    void updatehost();
    void updatedev();

inline   T &operator[](const GSIZET i) {
    #if defined(_G_BOUNDS_CHK)
      const char serr[] = "GTVector<T>::operator[]: ";
      if ( !bconstdata_ && i+gindex_.beg() > gindex_.end() ) {
        std::cout << serr << "Access error: " << i << std::endl;
        while(1){}; exit(1);
      }
    #endif
      return bconstdata_ ? data_[0] : data_[i+gindex_.beg()];
    };

inline    T operator[](const GSIZET i) const {
    #if defined(_G_BOUNDS_CHK)
      const char serr[] = "GTVector<T>::operator[] const: ";
      if ( !bconstdata_  && i+gindex_.beg() > gindex_.end() ) {
        std::cout << serr << "Access error: " << i << std::endl;
        while(1){}; exit(1);
      }
    #endif
      T ret = bconstdata_ ? data_[0] : data_[i+gindex_.beg()];
      return ret;
    };


    #pragma acc routine vector 
    GTVector<T>       &operator+(GTVector<T> &b);
    #pragma acc routine vector 
    GTVector<T>       &operator-(GTVector<T> &b);

    #pragma acc routine vector
    void               operator+=(T b);
    #pragma acc routine vector
    void               operator-=(T b);
    #pragma acc routine vector
    void               operator*=(T b);
    #pragma acc routine vector
    GTVector<T>       &operator*(T b);
  
    #pragma acc routine vector 
    void               operator+=(GTVector<T> &b);
    #pragma acc routine vector
    void               operator-=(GTVector<T> &b);
    #pragma acc routine vector
    void               operator*=(GTVector<T> &b);


    #pragma acc routine vector
    T max();
    #pragma acc routine vector
    T maxn(GSIZET n);
    #pragma acc routine vector
    T min();
    #pragma acc routine vector
    T minn(GSIZET n);
    #pragma acc routine vector
    T sum();
    #pragma acc routine vector
    T L1norm();
    #pragma acc routine vector
    T L2norm();
    #pragma acc routine vector
    GBOOL contains(T val);     // Buffer contains val?
    #pragma acc routine vector
    GSIZET contains(T val, GSIZET *&iwhere, GSIZET &nw);
    #pragma acc routine vector
    GBOOL containsn(T val, GSIZET n); // check n elements
    #pragma acc routine vector
    GBOOL  contains(T val, GSIZET  &index);
    #pragma acc routine vector
    GBOOL  containsn(T val, GSIZET n, GSIZET  &index);
    #pragma acc routine vector
    GBOOL  contains_floor(T val, GSIZET &index, T floor, GSIZET istart=0);
    #pragma acc routine vector
    GBOOL  contains_ceil(T val, GSIZET  &index, T ceil  , GSIZET istart=0);
    #pragma acc routine vector
    GSIZET distinctrng(GSIZET istart, GSIZET n, GSIZET is,  T *&vals, GSIZET *&index, GSIZET  &n_distinct);
    #pragma acc routine vector
    GSIZET distinctrng(GSIZET istart, GSIZET n, GSIZET is,  GSIZET *&index, GSIZET  &n_distinct);
    #pragma acc routine vector
    GSIZET distinctrng_floor(GSIZET istart, GSIZET n, GSIZET is, T *&vals, GSIZET *&index, GSIZET  &n_distinct, T floor);
    #pragma acc routine vector
    GSIZET distinctrng_floor(GSIZET istart, GSIZET n, GSIZET is, GSIZET *&index, GSIZET  &n_distinct, T floor);
    #pragma acc routine vector
    GSIZET distinct(GSIZET *&index, GSIZET  &n_distinct);
    #pragma acc routine vector
    GSIZET distinct_floor(GSIZET *&index, GSIZET  &n_distinct, T floor);
    #pragma acc routine vector
    GSIZET multiplicity(T val);
    #pragma acc routine vector
    GSIZET multiplicity(T val, GSIZET *&index, GSIZET &n);
    #pragma acc routine vector
    GSIZET multiplicity_floor(T val, T floor);
    #pragma acc routine vector
    GSIZET multiplicity_ceil(T val, T ceil);
    #pragma acc routine vector
    GSIZET multiplicity_floor(T val, GSIZET *&index, GSIZET &n, T floor);

    #pragma acc routine vector
    void   sortincreasing();
    #pragma acc routine vector
    void   sortincreasing(GTVector<GSIZET> &);
    #pragma acc routine vector
    void   sortdecreasing();
    #pragma acc routine vector
    void   sortdecreasing(GTVector<GSIZET> &);


  private:

    GIndex gindex_; // gen. index object
    T     *data_;
    GSIZET n_;
    GBOOL  bdatalocal_; // tells us that data_ is owned by caller
    GBOOL  bconstdata_; // says data is const, and access is to 1st element only


  #pragma acc routine vector 
  GTVector<T>& add_impl_(GTVector<T> &b, std::true_type);
  #pragma acc routine vector 
  GTVector<T>& add_impl_(GTVector<T> &b, std::false_type);
  
  #pragma acc routine vector 
  GTVector<T>& sub_impl_(GTVector<T> &b, std::true_type);
  #pragma acc routine vector 
  GTVector<T>& sub_impl_(GTVector<T> &b, std::false_type);

};


//
// Writing to ostream doesn't need access to internal data so
// don't bother making it a friend class
//
template<typename T>
std::ostream &operator<<(std::ostream &os, GTVector<T> &obj) {
  for ( GLONG j=0; j < obj.size(); j++ ) {
    os << obj[j] << " ";
  }
  return os;
};

typedef GTVector<GFTYPE>   GVector;
typedef GTVector<GINT>     GIVector;
typedef GTVector<GINT>     GIBuffer;
typedef GTVector<GNODEID>  GNIDBuffer;
typedef GTVector<GSIZET>   GSZBuffer;
typedef GTVector<GFLOAT>   GFVector;
typedef GTVector<GDOUBLE>  GDVector;
typedef GTVector<GDOUBLE>  GGVector;
typedef GTVector<GQUAD>    GQVector;

#include "gtvector.ipp"

#endif
