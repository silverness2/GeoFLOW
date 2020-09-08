//==================================================================================
// Module       : gstress.hpp
// Date         : 09/05/20 (DLR)
// Description  : Represents the SEM discretization of the full viscous
//                stress-energy operator. The viscous stress in the 
//                momentum eqution is
//                    2 [ mu s_{ij}],j,
//                where
//                    s_{ij} = (Del_i u_j + Del_j u_i)/2
//                and the viscous stress-energy for the energy equation is
//                    2 [mu u_i s_{ij} ],j
//                where u_i is the velocity, and mu, the viscosity. Repeated
//                indices are summed here.  For the energy, this is a nonlinear 
//                operator, so should not derive from GLinOp. Operator requires 
//                that grid consist of elements of only one type.
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none
//==================================================================================


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Default constructor
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename TypePack>
GStressEnOp<TypePack>::GStressEnOp(Grid &grid)
:
bown_mu_      (TRUE),
bown_kappa_   (TRUE),
grid_         (&grid),
massop_       (&grid.massop()),
mu_           (NULLPTR)
{
  assert(grid_->ntype().multiplicity(0) == GE_MAX-1 
        && "Only a single element type allowed on grid");
  mu_     = new GTVector<Ftype>(1);
  *mu_    = 1.0;
  kappa_  = new GTVector<Ftype>(1);
  *kappa_ = 1.0;
} // end of constructor method (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename TypePack>
GStressEnOp<TypePack>::~GStressEnOp()
{
  if ( mu_    != NULLPTR && bown_mu_    ) delete mu_;
  if ( kappa_ != NULLPTR && bown_kappa_ ) delete kappa_;
} // end, destructor


//**********************************************************************************
//**********************************************************************************
// METHOD : set_mu
// DESC   : Set viscosity (for momentum). This may be a field if length>1, 
//          or a constant if length == 1.
// ARGS   : 
//          mu    : 'viscosity' parameter, global
//             
// RETURNS:  none
//**********************************************************************************
template<typename TypePack>
void GStressEnOp<TypePack>::set_mu(StateComp &mu)
{
  assert(mu.size() == 1 || mu.size() >= grid_->ndof()
       && "Viscosity parameter of invalid size");

  if ( mu_ != NULLPTR && bown_mu_ ) delete mu_;

  mu_ = &mu;
  bown_mu_ = FALSE;

} // end of method set_mu


//**********************************************************************************
//**********************************************************************************
// METHOD : set_kappa
// DESC   : Set kappa (for stress energy). This may be a field if length>1, 
//          or a constant if length == 1.
// ARGS   : 
//          kappa : 'kappa' parameter, global
//             
// RETURNS:  none
//**********************************************************************************
template<typename TypePack>
void GStressEnOp<TypePack>::set_kappa(StateComp &kappa)
{
  assert(kappa.size() == 1 || kappa.size() >= grid_->ndof()
       && "Energy dissipation parameter of invalid size");

  if ( kappa_ != NULLPTR && bown_kappa_ ) delete kappa_;

  kappa_ = &kappa;
  bown_kappa_ = FALSE;

} // end of method set_kappa


//**********************************************************************************
//**********************************************************************************
// METHOD : apply (1)
// DESC   : Compute application of this operator to input momentum vector:
//            so = [ mu (Del_i u_j + Del_j u_i)],j
//          Remember, normally, this discretization would be multiplied by
//          -1 to represent this operation. We do not apply this sign here.
// ARGS   : u   : input vector field
//          idir: which momentum component we're computing for
//          utmp: array of tmp arrays
//          so  : output (result) vector component, idir
//             
// RETURNS:  none
//**********************************************************************************
template<typename TypePack>
void GStressEnOp<TypePack>::apply(State &u, GINT idir, State &utmp, StateComp &so) 
{

  assert( utmp.size() >= 3
       && "Insufficient temp space specified");

  GINT       nxy = grid_->gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;
  GElemList *gelems=&grid_->elems();

  assert( idir > 0 && idir <= nxy );


  // so = D^{T,j} [mu [D_i u_j + Dj u_i)]:

  // Do D^{T,j} [mu(D_j u_i) ] terms:
  grid_->deriv(*u[idir-1], 1, FALSE, *utmp[0], *utmp[1]);
  if ( mu_->size() > 1 ) { // point-multiply by mu:
    utmp[2]->pointProd(*mu_);
  }
  else {
    *utmp[2] *= (*mu_)[0];
  }
  grid_->deriv(*utmp[1]  , 1, TRUE , *utmp[0], so);
  for ( auto j=1; j<nxy; j++ ) { 
     grid_->deriv(*u[idir-1], j+1, FALSE, *utmp[0], *utmp[1]);
     // Point-multiply by mu before taking 'divergence':
     if ( mu_->size() > 1 ) {
       utmp[1]->pointProd(*mu_);
     }
     else {
       *utmp[1] *= (*mu_)[0];
     }
     grid_->deriv(*utmp[1]  , j+1, TRUE , *utmp[0], *utmp[2]);
     so += *utmp[2];
  }

  // Do D^{T,j} [mu(D_i u_j) ] terms:
  for ( auto j=0; j<nxy; j++ ) { 
     grid_->deriv(*u[j], idir, FALSE, *utmp[0], *utmp[1]);
     // Point-multiply by mu before taking 'divergence':
     if ( mu_->size() > 1 ) {
       utmp[1]->pointProd(*mu_);
     }
     else {
       *utmp[1] *= (*mu_)[0];
     }
     grid_->deriv(*utmp[1]  , j+1, TRUE , *utmp[0], *utmp[2]);
     so += *utmp[2];
  }

  // Point-multiply by mass & Jacobian:
  so *= *(massop_->data());

} // end of method apply (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : apply (2)
// DESC   : Compute application of this operator to input energy:
//            eo = [ mu u_i (Del_i u_j + Del_j u_i)],j
//          Remember, normally, this discretization would be multiplied by
//          -1 to represent this operation. We do not apply this sign here.
// ARGS   : u   : input vector field
//          utmp: array of tmp arrays
//          eo  : output (result) vector component, idir
//             
// RETURNS:  none
//**********************************************************************************
template<typename TypePack>
void GStressEnOp<TypePack>::apply(State &u, State &utmp, StateComp &eo) 
{

  assert( utmp.size() >= 3
       && "Insufficient temp space specified");

  GINT       nxy = grid_->gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;
  GElemList *gelems=&grid_->elems();


  // eo = D^{T,j} [mu u^i [D_i u_j + Dj u_i)]:

  // Do D^{T,j} [ mu u^i (D_j u_i) ] terms:
  eo = 0.0;
  for ( auto j=0; j<nxy; j++ ) { 
    *utmp[1] = 0.0;
    for ( auto i=0; i<nxy; i++ ) {
       grid_->deriv(*u[i], j+1, FALSE, *utmp[0], *utmp[2]);
       utmp[2]->pointProd(*u[i]);
       *utmp[1] += *utmp[1];
    }
    // Point-multiply by kappa before taking 'divergence':
    if ( kappa_->size() > 1 ) {
      utmp[1]->pointProd(*kappa_);
    }
    else {
      *utmp[1] *= (*kappa_)[0];
    }
    grid_->deriv(*utmp[1]  , j+1, TRUE , *utmp[0], *utmp[2]);
    eo += *utmp[2];
  }


  // Do D^{T,j} [ mu u^i (D_i u_j) ] terms:
  for ( auto j=0; j<nxy; j++ ) { 
    *utmp[1] = 0.0;
    for ( auto i=0; i<nxy; i++ ) {
       grid_->deriv(*u[j], i+1, FALSE, *utmp[0], *utmp[2]);
       utmp[2]->pointProd(*u[i]);
       *utmp[1] += *utmp[1];
    }
    // Point-multiply by kappa before taking 'divergence':
    if ( kappa_->size() > 1 ) {
      utmp[1]->pointProd(*kappa_);
    }
    else {
      *utmp[1] *= (*kappa_)[0];
    }
    grid_->deriv(*utmp[1]  , j+1, TRUE , *utmp[0], *utmp[2]);
    eo += *utmp[2];
  }

  // Point-multiply by mass & Jacobian:
  eo *= *(massop_->data());

} // end of method apply (2)

