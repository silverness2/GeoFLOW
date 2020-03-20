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
template<typename Types>
GCG<Types>::GCG(Traits& traits, Grid& grid, ConnectivityOp& ggfx, State& tmppack)
: SolverBase(traits, grid, ggfx, tmppack),
comm_        (ggfx.getComm()),
bInit_                (FALSE),
bbv_                  (FALSE),
nprocs_                   (1),
irank_                    (0),
precond_            (NULLPTR)
{
  irank_   = GComm::WorldRank(comm_);
  nprocs_ = GComm::WorldSize(comm_);

} // end of constructor (2) method


//************************************************************************************
//************************************************************************************
// METHOD : Destructor
// DESC   : 
// ARGS   : none.
// RETURNS: none.
//************************************************************************************
template<typename Types>
GCG<Types>::~GCG()
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
template<typename Types>
GCG<Types>::GCG(const GCG<Types> &a)
{

} // end of copy constructor method


//************************************************************************************
//************************************************************************************
// METHOD : Assignment operatior
// DESC   : 
// ARGS   : GCG
// RETURNS: GCG
//************************************************************************************
template<typename Types>
GCG<Types>& GCG<Types>::operator=(const GCG<Types> &a)
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
template<typename Types>
void GCG<Types>::init()
{
  if ( bInit_ ) return;

  residuals_.resize(this->traits_.maxit);
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
template<typename Types>
GINT GCG<Types>::solve_impl(Operator& A, const StateComp& b, StateComp& x)
{
  GINT       iret=GCGERR_NONE;
  GFTYPE     alpha, beta, residual, rho, rhom;
  StateComp *q, *r, *w, *z;
  State      tmp(this->tmp_->size()-4);
  StateComp *mask = &this->grid_->get_mask();

  assert(this->tmp_->size() > 5);
  init();

  // Set some pointers:
  tmp.resize(this->tmp_->size()-4);
  q  = (*this->tmp_)[1];
  r  = (*this->tmp_)[3];
  w  = (*this->tmp_)[0];
  z  = (*this->tmp_)[2];
  for ( auto j=0; j<this->tmp_->size()-4; j++ ) {
    tmp[j] = (*this->tmp_)[j+4];
  }
  residuals_ = 0.0;

 // Initialize CG loop, and enter loop:
 *r = b;
//if ( bbv_ ) x.pointProd(*mask);
  A.opVec_prod(x, tmp, *p);             // Ax
 *r -= (*p);                            // r = b - Ax, initial residual
cout << "GCG::solve: r=" << *r << endl;
  this->ggfx_->doOp(*r, GGFX_OP_SUM);   // DSS r
  if ( bbv_ ) r->pointProd(*mask);      // Mask DSS r
  if ( precond_ != NULLPTR ) {          // solve P z = r for z
    iret = precond_->solve(*r, *z);  
    if ( iret >  0 ) iret = GCGERR_PRECOND; 
  }
  else {
    *z = *r;                            // use identity preconditioner
  }
  iter_ = 0; residual = 1.0;

  while ( iret == GCERR_NONE 
       && iter_ < this->traits_.maxit 
       && residual > this->traits_.tol ) {

    if ( precond_ != NULLPTR ) {        // solve P z = r for z
      iret = precond_->solve(*r, *z);   // apply preconditioner
      if ( iret >  0 ) {
        iret = GCGERR_PRECOND;
        break;
      }
    }
    else {
      *z = *r;                          // use identity preconditioner
    }
//cout << "GCG::solve: iter=" << iter_ << " zk=" << *z << endl;
    rho = r->gdot(*z,comm_);
//cout << "GCG::solve: iter=" << iter_ << " rho=" << rho << " rhom=" << rhom << endl;
    if ( iter_  == 0 ) {                // find p
     *p = *z;
    }
    else {
      beta = rho / rhom;
//cout << "GCG::solve: iter=" << iter_ << " beta =" << beta << endl;
      GMTK::saxpby(*p, beta, *z, 1.0);
    }

    A.opVec_prod(*p, tmp, *q);          // q = A p
    this->ggfx_->doOp(*q, GGFX_OP_SUM); // DSS q
//cout << "GCG::solve: qk=" << *q << endl;
    alpha = rho / (p->gdot(*q,comm_));
//cout << "GCG::solve: iter=" << iter_ << " alpha=" << alpha << endl;
//  if ( bbv_ ) p->pointProd(*mask); 
    GMTK::saxpby( x, 1.0, *p, alpha);   // x = x + alpha p
    GMTK::saxpby(*r, 1.0, *q,-alpha);   // r = r - alpha q
//  if ( bbv_ ) r->pointProd(*mask);

    rhom = rho;

    residual = compute_norm(*r, tmp); // find norm of residual
//cout << "GCG::solve: residual=" << residual << " tol=" << this->traits_.tol << endl;
    residuals_[iter_] = residual;
    iter_++;

  } // end, CG loop

  if ( bbv_ ) x.pointProd(*mask);

  cout << "GCG::solve: iter_     =" << iter_ << endl;
  cout << "GCG::solve: rersiduals=" << residuals_ << endl;
    
  if ( iret == GCGERR_NONE 
    && iter_ >= this->traits_.maxit 
    && residual > this->traits_.tol ) iret = GCGERR_NOCONVERGE;

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
template<typename Types>
GINT GCG<Types>::solve_impl(Operator& A, const StateComp& b, const StateComp& xb, StateComp& x)
{
  GINT iret;

  bbv_ = TRUE;

  // Add in boundary solution, so 
  // that it will form RHS in homogeneous solve:

  if ( bbv_ ) x.pointProd(this->grid_->get_mask()); // apply mask
  x += xb; 

  // Compute homogeneous solution:
  iret = solve_impl(A, b, x);
  assert(iret == GCGERR_NONE);

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
template<typename Types>
GFTYPE GCG<Types>::compute_norm(StateComp& x, State& tmp)
{
  GFTYPE lret[2], gret[2];

  switch (this->traits_.normtype) {
    case LinSolverBase<Types>::GCG_NORM_INF:
      lret[0] = x.infnorm();
      GComm::Allreduce(lret, gret, 1, T2GCDatatype<GFTYPE>() , GC_OP_MAX, comm_);
      break;
    case LinSolverBase<Types>::GCG_NORM_EUC:
     *tmp[0]  = x; tmp[0]->rpow(2);
      lret[0]  = tmp[0]->sum();
      lret[1]  = x.size();
      GComm::Allreduce(lret, gret, 2, T2GCDatatype<GFTYPE>() , GC_OP_SUM, comm_);
      gret[0]  = sqrt(gret[0]/gret[1]);
      break;
    case LinSolverBase<Types>::GCG_NORM_L2:
     *tmp[0]  = x; tmp[0]->rpow(2);
      gret[0]  = this->grid_->integrate(*tmp[0],*tmp[1]);
      gret[0] /= this->grid_->volume();
      break;
    case LinSolverBase<Types>::GCG_NORM_L1:
     *tmp[0]  = x; tmp[0]->abs();
      gret[0]  = this->grid_->integrate(*tmp[0],*tmp[1]);
      gret[0] /= this->grid_->volume();
      break;
    default:
      assert(FALSE);
  }

  return gret[0];

} // end, compute_norm


