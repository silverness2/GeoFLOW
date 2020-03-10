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

    

} // end of method solve_impl (1)


//************************************************************************************
//************************************************************************************
// METHOD : solve_impl (2)
// DESC   : Solve implementation, with boundary solution
// ARGS   : A    : linear operator to invert
//          b    : right-hand side vector
//          xb   : boundary solution
//          x0   : homogeneous solution, returned
// ARGS   : 
// RETURNS: integer error code
//************************************************************************************
template<typename TPack>
GINT GCG<T>::solve_impl(Operator& A, const StateComp& b, StateComp& xb, StateComp& x0)
{
  GINT      n=0;
  GFTYPE    alpha, rho, rhom;
  StateComp *p, *q, *r, *rm, *zm;

  p  = &(this->tmp_)[0];
  q  = &(this->tmp_)[1];
  zm = &(this->tmp_)[2];
  rm = &(this->tmp_)[3];
  r  = &(this->tmp_)[4];

  A.opVec_prod(x, utmp, *p);
  *r = b - (*p); // initial residual
  while ( n<traits_.maxit && residual > traits_.tol ) {
    if ( precond_ != NULLPTR ) {
      precond_->solve(*r, *z); // apply preconditioner
    }
    else {
      *z = *r;                 // use identity preconditioner
    }
  } // end, CG loop

} // end of method solve_impl (2)


