//==================================================================================
// Module       : gboyd_filter.hpp
// Date         : 9/14/20 (DLR)
// Description  : Computes the Boyd filter to diminish aliasing errors.
//                Taken from Giraldo & Rosemont 2004, MWR:132 133:
//                    u <-- F u
//                where
//                    F = L Lambda L^-1; s.t.
//                and 
//                    Lambda = 1 if i< ifilter
//                             mu [(i-ifilter)/(N - ifilter)]^2, i>= ifilter.
//                L is the Legendre transform matrix:
//                    L = | P_0(xi0), P_1(xi0) ... P_i(xi0)-P_{i-2)(xi0) ... |
//                        | P_0(xi1), P_1(xi1) ... P_i(xi1)-P_{i-2)(xi1) ... |
//                        |  ...                                             |.
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : FilterBase
//==================================================================================

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Default constructor
// ARGS   : grid    : Grid object
//          ifilter : starting mode for filtering
//          mufilter: filter factor (by which to truncatte)
// RETURNS: none
//**********************************************************************************
template<typename TypePack>
GBoydFilter<TypePack>::GBoydFilter(Traits &traits, Grid &grid)
:
bInit_               (FALSE),
ifilter_    (traits.ifilter),
mufilter_  (traits.mufilter),
grid_                (&grid)
{
  assert(grid_->ntype().multiplicity(0) == GE_MAX-1 
        && "Only a single element type allowed on grid");
} // end of constructor method (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename TypePack>
GBoydFilter<TypePack>::~GBoydFilter()
{
} // end, destructor


//**********************************************************************************
//**********************************************************************************
// METHOD : apply_impl
// DESC   : Compute application of this filter to input vector
//           
// ARGS   : t   : Time
//          u   : input vector field
//          utmp: array of tmp arrays; not used here
//          uo  : output (result) vector
//             
// RETURNS:  none
//**********************************************************************************
template<typename TypePack>
void GBoydFilter<TypePack>::apply_impl(const Time &t, StateComp &u, State &utmp, StateComp &uo) 
{

  GSIZET           ibeg, iend; // beg, end indices in global array
  GTVector<Ftype>  tmp;
  GTMatrix<Ftype> *F[GDIM];
  GElemList       *gelems=&grid_->elems();

  if ( !bInit_ ) init();

  for ( auto e=0; e<gelems->size(); e++ ) {
    ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
    u .range(ibeg, iend); // restrict global vecs to local range
    uo.range(ibeg, iend); 
    F [0] = (*gelems)[e]->gbasis(0)->getFilterMat();
    F [1] = (*gelems)[e]->gbasis(1)->getFilterMat(TRUE);
#if defined(_G_IS2D)
    GMTK::D2_X_D1(*F[0], *F[1], u, tmp, uo);
#elif defined(_G_IS3D)
    F [2] = (*gelems)[e]->gbasis(2)->getFilterMat(TRUE);
    GMTK::D3_X_D2_X_D1(*F[0], *F[1], *F[2],  u, tmp, uo);
#endif
  }
  u .range_reset(); 
  uo.range_reset(); 

} // end of method apply_impl


//**********************************************************************************
//**********************************************************************************
// METHOD : apply_impl
// DESC   : In-place application of this filter to input vector
//           
// ARGS   : t   : Time
//          u   : input vector field
//          utmp: array of tmp arrays
//             
// RETURNS:  none
//**********************************************************************************
template<typename TypePack>
void GBoydFilter<TypePack>::apply_impl(const Time &t, StateComp &u, State &utmp) 
{

  assert( utmp.size() >= 1
       && "Insufficient temp space provided");

  apply_impl(t, u, utmp, *utmp[0]); 
  u = *utmp[0];

} // end of method apply_impl


//**********************************************************************************
//**********************************************************************************
// METHOD : init
// DESC   : Initilize operators
// ARGS   : none.
// RETURNS:  none
//**********************************************************************************
template<typename TypePack>
void GBoydFilter<TypePack>::init()
{

  GINT             nnodes;
  GTMatrix<Ftype> *F, *FT, *iL, *L, Lambda;
  GTMatrix<Ftype> tmp;
  GElemList       *gelems=&grid_->elems();

  // Build the filter matrix, F, and store within 
  // each basis object for later use:
  //   F = L Lambda L^-1; s.t.
  // u <-- F u
  // where
  //   Lambda = 1 if i< ifilter
  //            mu [(i-ifilter)/(N - ifilter)]^2, i>= ifilter
  for ( auto e=0; e<gelems->size(); e++ ) {
    for ( auto k=0; k<GDIM; k++ ) {
      F  = (*gelems)[e]->gbasis(k)->getFilterMat(); // storage for filter
      FT = (*gelems)[e]->gbasis(k)->getFilterMat(TRUE); // transpose 
      nnodes = (*gelems)[e]->gbasis(k)->getOrder()+1;
      Lambda.resize(nnodes,nnodes); Lambda = 0.0;
      tmp   .resize(nnodes,nnodes);
      L      = (*gelems)[e]->gbasis(k)->getLegTransform();
      iL     = (*gelems)[e]->gbasis(k)->getiLegTransform();
      for ( auto i=0; i<nnodes; i++ ) { // build weight matix, Lambda
        Lambda(i,i) = 1.0;
        if ( i >= ifilter_ ) {
          Lambda(i,i) = mufilter_
                      * pow( ( (Ftype)(i-ifilter_) / ( (Ftype)(nnodes-ifilter_) ) ), 2);
        }
      } // end, node/mode loop 
      tmp = Lambda * (*iL);
     *F   = (*L) * tmp;
      F   ->transpose(*FT);
    } // end, k-loop
  } // end, element loop

  bInit_ = TRUE;

} // end of method init


