//==================================================================================
// Module       : gtpoint
// Date         : 7/1/18 (DLR)
// Description  : Encapsulates the access methods and data associated with
//                defining template 'point' object.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#include "gtypes.h"
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <typeinfo>
#include "gtvector.hpp"


template <typename T> class GTPoint;
template<typename T> std::ostream &operator<<(std::ostream &, GTPoint<T> &);

#if !defined(GTPOINT_HPP)
#define GTPOINT_HPP

template<typename T> class GTPoint
{
private:
  GINT           gdim_;
  GTVector<GINT> nexp_;
  T              eps_;
  GTVector<T>    lr_;
  GTVector<T*>   px_;

  template<typename U=T> typename std::enable_if<std::is_floating_point<U>::value, void>::type
  decompose()
  { for ( GINT i=0; i<gdim_; i++ ) {
    lr_[i] = roundl(frexp((*px_[i]), &nexp_[i])* 1.0e8); } }

public:
                  GTPoint(); 
                  GTPoint(GINT dim); 
                  GTPoint(GINT idim, T eps); 
                  GTPoint(const GTPoint<T> &e); 
                 ~GTPoint();

  T x1;
  T x2;
  T x3;
  T x4;

         void      resize(GINT dim);
         void      setBracket(T eps);
         T         getBracket();
  inline GTVector<T*> &data() 
  { return px_;}

  inline GINT     size() 
  { return gdim_;}
  inline GINT     dim() 
  { return gdim_;}


  template<typename U=T> typename std::enable_if<std::is_floating_point<U>::value, GBOOL>::type
  operator==(const GTPoint<T> &pp) // Add 'fuzziness' to equality check
  { GINT  i, n; GBOOL b; T r, del; GTPoint pq=pp;    
    // p.x_i = g * 2^n = g * 10^m, where m=n log_10(2):
    decompose();
    for ( i=0,b=TRUE;i<gdim_;i++) { r=frexp(pq[i],&n); ; // scale bracket to point
      b = b && ( nexp_[i]==n && lr_[i]==lround(r*1.0e8) ); }      // do check for equality
    return b; }

  template<typename U=T> typename std::enable_if<!std::is_floating_point<U>::value, GBOOL>::type
  operator==(const GTPoint<T> &pp) // Add 'fuzziness' to equality check
  { GBOOL b; GTPoint<T> pq = pp;
    for ( GINT i=0,b=TRUE;i<gdim_;i++) b = b && ( pq[i] == *px_[i] );
    return b; }

  template<typename U=T> typename std::enable_if<std::is_floating_point<U>::value, T>::type
  norm()  // Euclidean norm
  { T xn=0; for (GINT i=0; i<gdim_; i++ ) xn += ((*px_[i])*(*px_[i])); return sqrt(xn); }

  inline GBOOL    operator!=(const GTPoint<T> &p) 
  { return !this->operator==(p); }

  inline void  operator=(const GTPoint<T> &p)
  { eps_ = p.eps_; x1 = p.x1; x2 = p.x2; x3 = p.x3; x4 = p.x4;}

  inline GTPoint<T> &operator=(T f)
  { x1 = f;  x2 = f;  x3 = f;  x4 = f; return *this;}

  inline void  operator-(T f)
  { GTPoint<T> *ret = new GTPoint<T>(gdim_);
    for ( GINT j=0; j<gdim_; j++ ) (*ret)[j] = *px_[j] - f; return *ret; }

  inline GTPoint<T> &operator-(GTPoint<T> &a)
  { GTPoint<T> *ret=new GTPoint<T>(gdim_);
    for ( GINT j=0; j<gdim_; j++ ) (*ret)[j] = *px_[j] - a[j];  return *ret; }

  inline void  operator-=(T a)
  { for ( GINT j=0; j<gdim_; j++ ) *px_[j] -= a; }

  inline void  operator-=(GTPoint<T> &a)
  { for ( GINT j=0; j<gdim_; j++ ) *px_[j] -= a[j]; }

  inline void  operator+(T f)
  { GTPoint<T> *ret = new GTPoint<T>(gdim_);
    for ( GINT j=0; j<gdim_; j++ ) (*ret)[j] = *px_[j] + f; return *ret; }

  inline GTPoint<T> &operator+(GTPoint<T> &a)
  { GTPoint<T> *ret=new GTPoint<T>(gdim_);
    for ( GINT j=0; j<gdim_; j++ ) (*ret)[j] = *px_[j] + a[j];  return *ret; }

  inline void  operator+=(T a)
  { for ( GINT j=0; j<gdim_; j++ ) *px_[j] += a; }

  inline void  operator+=(GTPoint<T> &a)
  { for ( GINT j=0; j<gdim_; j++ ) *px_[j] += a[j]; }

  inline GTPoint<T> &operator*(T a)
  { GTPoint<T> *ret = new GTPoint<T>(gdim_);
    for ( GINT j=0; j<gdim_; j++ ) (*ret)[j] = *px_[j] * a; return *ret; }

  inline void  operator*=(T a)
  { for ( GINT j=0; j<gdim_; j++ ) *px_[j] *= a; }

  inline void  operator/=(T a)
  { for ( GINT j=0; j<gdim_; j++ ) *px_[j] /= a; }

  inline T &operator()(const GINT i)
  { 
#if defined(_G_BOUNDS_CHK)
    if ( i<0 || i>=gdim_ ) {
      while(1);  
      std::cout << "GTPoint::(): access error; bad index: " << i << std::endl; exit(1); }
#endif
    return *px_[i]; }

  inline T &operator[](const GINT i)
  { 
#if defined(_G_BOUNDS_CHK)
    if ( i<0 || i>=gdim_ ) { 
      while(1);
      std::cout << "GTPoint::[]: access error; bad index: " << i << std::endl; exit(1); }
#endif
    return *px_[i]; }

  template<typename U=T> typename std::enable_if<std::is_floating_point<U>::value, void >::type
  truncate()
  { GINT  i, n; T del, is;
    // p.x_i = g * 2^n = g * 10^m, where m=n log_10(2):
    for ( i=0;i<gdim_;i++ ) { frexp(*px_[i],&n); del=eps_*pow(10.0,0.30103*n); // scale bracket to point
      is = *px_[i]>0?1.0:-1.0; *px_[i] = *px_[i]!=0.0 ? is*(fabs(*px_[i])-del): 0.0;} }

  inline T mag() { T v=0.0; for ( GINT j=0; j<gdim_; j++ ) v += pow(*px_[j],2.0); return sqrt(v); }

  friend std::ostream &operator<<(std::ostream &str, GTPoint<T> &obj) {
    str << "("  << *obj.px_[0];
    for ( GINT i=1; i<obj.gdim_; i++ ) {
    str << ", " << *obj.px_[i];} str << ")";
    return str;
  } // end, operator<<


};

typedef GTPoint<GFTYPE>  GFPoint;
typedef GTPoint<GINT>    GIPoint;

#include "gtpoint.ipp"

#endif