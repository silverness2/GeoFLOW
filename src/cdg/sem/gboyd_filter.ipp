//==================================================================================
// Module       : gboyd_filter.hpp
// Date         : 9/14/20 (DLR)
// Description  : Computes the Boyd filter to diminish aliasing errors.
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
GBoydFilter<TypePack>::GBoydFilter(Grid &grid, GINT ifilter, Ftype mufilter)
:
bInit_         (FALSE),
ifilter_       (ifilter),
mufilter_      (mufilter),
grid_          (&grid)
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
// METHOD : apply
// DESC   : Compute application of this filter to input vector
//           
// ARGS   : t   : Time
//          u   : input vector field
//          utmp: array of tmp arrays
//          uo  : output (result) vector
//             
// RETURNS:  none
//**********************************************************************************
template<typename TypePack>
void GBoydFilter<TypePack>::apply(Time t, StateComp &u, State &utmp, StateComp &uo) 
{

  assert( utmp.size() >= 2
       && "Insufficient temp space specified");

  GSIZET           ibeg, iend; // beg, end indices for global array
  GTMatrix<Ftype> *lambda(GDIM), *L(GDIM), *iL(GDIM);
  GElemList       *gelems=&grid_->elems();

  if ( !bInit ) init();

  for ( auto e=0; e<gelems->size(); e++ ) {
    ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
    u.range(ibeg, iend); // restrict global vecs to local range
    for ( auto k=0; k<GDIM; k++ ) {
      N     [k]= (*gelems)[e]->size(k);
      lambda[k] = (*gelems)[e]->gbasis(k)->getLegTransWeightMat();
      L     [k] = (*gelems)[e]->gbasis(k)->getLegTransform();
      iL    [k] = (*gelems)[e]->gbasis(k)->getiLegTransform();
    }
#if defined(_G_IS2D)
    GMTK::D2_X_I1(*Di, u, N[0], N[1], du);
#elif defined(_G_IS3D)
#endif
    u.range_reset(); 
  }

} // end of method apply


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
  GTMatrix<Ftype> *lambda;
  GElemList *gelems=&grid_->elems();

  // Build the weighting matrix, and store within 
  // each basis object:
  for ( auto e=0; e<gelems->size(); e++ ) {
    for ( auto k=0; k<GDIM; k++ ) {
      lambda = (*gelems)[e]->gbasis(k)->getLegTransWeightMat();
      nnodes = (*gelems)[e]->gbasis(k)->getOrder()+1;
     *lambda = 0.0;
      for ( auto i=0; i<nnodes; i++ ) {
        (*lambda)(i,i) = 1.0;
        if ( i >= ifilter_ ) {
          (*lambda)(i,i) = mufilter_
                         * pow((i-ifilter_)/(nnodes-ifilter_),2);
        }
      } // end, node/mode loop 
    } // end, k-loop
  } // end, element loop

  bInit_ = TRUE;

} // end of method init


