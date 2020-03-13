//==================================================================================
// Module       : gcg.hpp
// Date         : 3/7/20 (DLR)
// Description  : Encapsulates the methods and data associated with
//                a conjugate Gradient (CG) solver
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#include <cmath>
#include <type_traits>
#include <cassert>
#include <limits>
#include "ggfx.hpp"
#include "gcomm.hpp"
#include "tbox/error_handler.hpp"

using namespace std;

//************************************************************************************
//************************************************************************************
// METHOD : Constructor
// DESC   : 
// ARGS   : GD_COMM object, defaults to GD_COMM_WORLD
// RETURNS: GCG
//************************************************************************************
template<typename TPack>
GCG<T>::GCG(Traits& traits, Grid& grid, ConnectivityOp& ggfx, State& tmppack)
: LinSolverBase(Traits& traits, Grid& grid, ConnectivityOp& ggfx, State& tmppack)
comm_        (ggfx.getComm()),
bInit_                (FALSE),
nprocs_                   (1),
irank_                    (0),
precond_            (NULLPTR)
{
  rank_   = GComm::WorldRank(comm_);
  nprocs_ = GComm::WorldSize(comm_);

} // end of constructor (2) method


//************************************************************************************
//************************************************************************************
// METHOD : Destructor
// DESC   : 
// ARGS   : none.
// RETURNS: none.
//************************************************************************************
template<typename TPack>
GCG<T>::~GCG()
{
}


//************************************************************************************
//************************************************************************************
// METHOD : Copy constructor
// DESC   : 
// ARGS   : GCG
// RETURNS: none
//************************************************************************************
// Copy constructor method
template<typename TPack>
GCG<T>::GCG(const GCG<T> &a)
{

} // end of copy constructor method


//************************************************************************************
//************************************************************************************
// METHOD : Assignment operatior
// DESC   : 
// ARGS   : GCG
// RETURNS: GCG
//************************************************************************************
template<typename TPack>
GCG<T>  &GCG<T>::operator=(const GCG<T> &a)
{

  return *this;
 
} // end of = operator


//************************************************************************************
//************************************************************************************
// METHOD : init
// DESC   : Performs initialization of global GCG operator
// ARGS   : none.
// RETURNS: none.
//************************************************************************************
template<typename TPack>
void GCG<T>::init(GNIDBuffer &glob_index)
{
  residuals_.resize(traits_.maxit);
  bInit_ = TRUE;

} // end of method init


//************************************************************************************
//************************************************************************************
// METHOD : solve_impl (1)
// DESC   : Solve implementation
// ARGS   : A    : linear operator to invert
//          b    : right-hand side vector
//          x    : solution, returned
// RETURNS: integer error code
//************************************************************************************
template<typename TPack>
GINT GCG<T>::solve_impl(Operator& A, const StateComp& b, StateComp& x)
{
  GINT       iret;
  GFTYPE     alpha, residual, rho, rhom;
  StateComp *p, *q, *r, *z;
  State      tmp(this->tmp_->size()-4);

  assert(this_->tmp_->size() > 5);

  iret = GCGERR_NONE;


  tmp(.resize(this->tmp_->size()-4);
  p  = (*this->tmp_)[0];
  q  = (*this->tmp_)[1];
  z  = (*this->tmp_)[2];
  r  = (*this->tmp_)[3];
  for ( auto j=0; j<this->tmp_->size()-4) {
    tmp[j] = (*this->tmp_)[j+4];
  }

 *r = b;
  A.opVec_prod(x, tmp, *p);
 *r -= (*p);                         // initial residual
  ggfx_->doOp(*r, GGFX_OP_SSUM);      // DSS r
  iter_ = 0; residual = 1.0;
  while ( iter_ < traits_.maxit && residual > traits_.tol ) {

    if ( precond_ != NULLPTR ) {      // solve Mz = r for z
      iret = precond_->solve(*r, *z); // apply preconditioner
      if ( iret >  0 ) {
        iret = GCGERR_PRECOND;
        break;
      }
    }
    else {
      *z = *r;                        // use identity preconditioner
    }
    rho = r->gdot(*z);
    if ( iter_  == 0 ) {              // find p
     *p = *z;
    }
    else {
      beta = rho / rhom;
      GMTK::saxpby(*p, beta, *z, 1.0);
    }

    A.opVec_prod(*p, tmp, *q);        // q = A p
    ggfx_->doOp(*q, GGFX_OP_SUM);     // DSS q
    alpha = rho / (p->gdot(*q));
    GMTK::saxpby(*x, 1.0, *p, alpha); // x = x + alpha p
    GMTK::saxpby(*r, 1.0, *q,-alpha); // r = r - alpha q

    rhom = rho;

    residual = compute_norm(*r, tmp); // find norm of residual

    residuals_[iter_] = residual;
    iter_++;

  } // end, CG loop
    
  if ( iret != GCGERR_NONE 
    && iter_ >= traits_.maxit 
    && residual > traits_.tol ) iret = GCGERR_NOCONVERGE;

  return iret;

} // end of method solve_impl (1)


//************************************************************************************
//************************************************************************************
// METHOD : solve_impl (2)
// DESC   : Solve implementation, with boundary solution
// ARGS   : A    : linear operator to invert
//          b    : right-hand side vector
//          xb   : boundary solution
//          x    : solution, x0+xb, returned, where x0=homogeneous solution
// ARGS   : 
// RETURNS: integer error code
//************************************************************************************
template<typename TPack>
GINT GCG<T>::solve_impl(Operator& A, const StateComp& b, StateComp& xb, StateComp& x)
{

  // Subtract out the boundary solution:
  x -= xb; 

  // Compute homogeneous solution:
  iret = solve_impl(A, b, x);
  assert(iret != GCGERR_NONE);

  // Add back boundary solution:
  x += xb;
 
  return iret;

} // end of method solve_impl (2)


//************************************************************************************
//************************************************************************************
// METHOD : compute_norm
// DESC   : Compute required norm of input vector, based
//          on traits
// ARGS   : x    : vector
//          tmp  : tmp vector
// RETURNS: required norm 
//************************************************************************************
template<typename TPack>
GINT GCG<T>::compute_norm(const StateComp& x, State& tmp)
{
  GFTYPE lret[2], grert[2];

  switch (traits_.normtype) {
    case GCG_NORM_INF:
      lret[0] = x.infnorm();
      GComm::Allreduce(lret, gret, 1, T2GCDatatype<GFTYPE>() , GC_OP_MAX, comm_);
      break;
    case GCG_NORM_EUCLID:
     *utmp[0]  = x; utmp[0]->pow(2);
      lret[0]  = utmp[0]->sum();
      lret[1]  = x.size();
      GComm::Allreduce(lret, gret, 2, T2GCDatatype<GFTYPE>() , GC_OP_SUM, comm_);
      gret[0]  = sqrt(gret[0]/gret[1]);
      break;
    case GCG_NORM_L2:
     *utmp[0]  = x; utmp[0]->pow(2);
      gret[0]  = grid_->integrate(*tmp[0],*tmp[1]);
      gret[0] /= this->grid_->volume();
      break;
    case GCG_NORM_L1:
     *utmp[0]  = x; utmp[0]->abs();
      gret[0]  = grid_->integrate(*tmp[0],*tmp[1]);
      gret[0] /= this->grid_->volume();
      break;
  }

  return gret[0];

} // end, compute_norm


