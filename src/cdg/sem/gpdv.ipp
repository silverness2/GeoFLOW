//==================================================================================
// Module       : gpdv.ipp
// Date         : 11/11/18 (DLR)
// Description  : Represents the SEM discretization of the 'pdV' operator:
//                p Div u. This is a nonlinear operator, so should not derive 
//                from GLinOp. This operator requires that grid consist of
//                elements of only one type.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved.
// Derived From : GLinOp
//==================================================================================

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Default constructor
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename TypePack>
GpdV<TypePack>::GpdV(Grid &grid)
:
grid_           (&grid),
massop_       (&grid.massop())
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
GpdV<TypePack>::~GpdV()
{
} // end, destructor


//**********************************************************************************
//**********************************************************************************
// METHOD : apply
// DESC   : Compute application of this operator to input vector:
//            po = p Div u
// ARGS   : p   : input p field
//          u   : input vector field
//          utmp: array of tmp arrays
//          po  : output (result) vector
//             
// RETURNS:  none
//**********************************************************************************
template<typename TypePack>
void GpdV<TypePack>::apply(StateComp &p, State &u, State &utmp, StateComp &po) 
{

  assert( utmp.size() >= 2
       && "Insufficient temp space specified");

  GElemList *gelems=&grid_->elems();

// Must compute:
//    po = p Div u 
//
  StateComp *Jac = &grid_->Jac();

  // Compute po += Gj (D^j u_j): 
  grid_->deriv(*u[0], 1, *utmp[0], po );
  for ( auto j=1; j<GDIM; j++ ) { 
     grid_->deriv(*u[j], TRUE, *utmp[0], *utmp[1] );
     po += *utmp[1];
  }

  // Point-multiply by p:
  if ( p.size() > 1 ) {
    po.pointProd(p);
  }
  else {
    po *= p[0];
  }

  // apply mass:
#if 1
  po *= *massop_->data();
  po *= *Jac;
#endif

} // end of method apply


