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
bInitialized_   (FALSE),
grid_           (&grid),
massop_       (&grid.massop())
{
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
  assert(bInitialized_ && "Operator not initialized");
    
  if ( grid_->itype(GE_REGULAR).size() > 0 ) {
    reg_prod(p, u, utmp, po);
  }
  else if ( grid_->itype(GE_DEFORMED).size() > 0 
         || grid_->itype(GE_2DEMBEDDED).size() > 0 ) {
    def_prod(p, u, utmp, po);
  }

} // end of method apply


//**********************************************************************************
//**********************************************************************************
// METHOD : def_prod
// DESC   : Compute this term for 2d & 3d GE_DEFORMED elements.
//          NOTE: must have 2 utmp vectors provided
// ARGS   : p   : input vector (global)
//          u   : input vector field (global)
//          utmp: array of tmp arrays. 2 arrays are required.
//          po  : output (result) field (global)
//             
// RETURNS:  none
//**********************************************************************************
template<typename TypePack>
void GpdV<TypePack>::def_prod(StateComp &p, State &u, State &utmp, StateComp &po) 
{
  assert( utmp.size() >= GDIM+1
       && "Insufficient temp space specified");

  GElemList *gelems=&grid_->elems();

// Must compute:
//    po = p Div u = p (G1 D1 + G2 D2 + G3 D3 ) u
//
// where
//    D1u = (I_X_I_X_Dx)u; D2u = (I_X_Dy_X_I)u; D3u = (Dz_X_I_X_I)u
// and Gj are the 'metric' terms computed in the element, dXi^j/dX^j
// that include the weights and Jacobian

  // Compute po += Gj (D^j u_j): 
  po = 0.0;
  for ( auto j=0; j<GDIM; j++ ) { 
    grid_->compute_grefdiv(u, etmp1_, FALSE, *utmp[0]); 
    utmp[0]->pointProd(*G_[j],*utmp[1]); // Gj * du^j. Mass included in Gj
    po += *utmp[1];
  }

  // Point-multiply by p:
  if ( p.size() > 1 ) {
    po.pointProd(p,*utmp[0]);
  }
  else {
    po *= p[0];
  }

} // end of method def_prod


//**********************************************************************************
//**********************************************************************************
// METHOD : reg_prod
// DESC   : Compute application of this operator to input vector for 
//          GE_REGULAR elements
//          NOTE: must have 1 utmp vectors provided
// ARGS   : p   : input vector (global)
//          u   : input vector field (global)
//          utmp: array of tmp arrays. One array required.
//          po  : output (result) field (global)
//             
// RETURNS:  none
//**********************************************************************************
template<typename TypePack>
void GpdV<TypePack>::reg_prod(StateComp &p, State &u, State &utmp, StateComp &po) 
{

  assert( utmp.size() >= GDIM+1
       && "Insufficient temp space specified");

  GElemList *gelems=&grid_->elems();

// Must compute:
//    po = p Div u = p (G1 D1 u1 + G2 D2 u2 + G3 D3 u3 )
//
// where
//    D1u = (I_X_I_X_Dx)u; D2u = (I_X_Dy_X_I)u; D3u = (Dz_X_I_X_I)u
// and Gj are the 'metric' terms computed in the element, dXi^j/dX^j
// that include the weights and Jacobian. For regular elements, 
// Gj are constant in each direction (but necessarily all equal), 
// as is the Jacobian. Weights are provided by mass matrix, which
// is set on construction, so that we don't have to re-compute it here:

  grid_->compute_grefdiv(u, etmp1_, FALSE, *utmp[0]); // utmp stores tensor-prod divergence: Dj u^j

  // Compute po += Gj (D^j u_j): 
  po = 0.0;
  for ( auto j=0; j<GDIM; j++ ) { 
    grid_->compute_grefdiv(u, etmp1_, FALSE, *utmp[0]); 
    *utmp[0] *= (*G_[j])[0]; // remember, mass not included in G here
    po += *utmp[0];
  }

  // Point-multiply by p:
  if ( p.size() > 1 ) {
    po.pointProd(p,*utmp[0]);
  }
  else {
    po *= p[0];
  }

  // apply mass:
  massop_->opVec_prod(*utmp[0], utmp, po); // last arg is unused

} // end of method reg_prod


//**********************************************************************************
//**********************************************************************************
// METHOD : init
// DESC   : Compute metric components for all elements. Split up this way
//          in case additional data must be initialized.
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename TypePack>
void GpdV<TypePack>::init()
{
  assert(grid_->ntype().multiplicity(0) == GE_MAX-1 
        && "Only a single element type allowed on grid");

  if ( grid_->itype(GE_2DEMBEDDED).size() > 0 
    || grid_->itype  (GE_DEFORMED).size() > 0 ) {
    def_init();
  }
  if ( grid_->itype(GE_REGULAR).size() > 0 ) {
    reg_init();
  }

  bInitialized_ = TRUE;

} // end, method init


//**********************************************************************************
//**********************************************************************************
// METHOD : def_init
// DESC   : Compute metric components if there are deformed elements.
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename TypePack>
void GpdV<TypePack>::def_init()
{

  if ( grid_->itype(GE_2DEMBEDDED).size() == 0 
    && grid_->itype  (GE_DEFORMED).size() == 0 ) return;


  GTVector<GSIZET>             N(GDIM);
  GTMatrix<StateComp>  *dXidX;    // element-based dXi/dX matrix
  GTVector<StateComp*>  W(GDIM);  // element-based weights
  StateComp            *Jac;      // element-based Jacobian
  GElemList                   *gelems = &grid_->elems();

  // Compute 'metric' components:
  // Gi = [dxi/dx, deta/dy, dzeta/dz]; 
  GSIZET nxy = grid_->itype(GE_2DEMBEDDED).size() > 0 ? GDIM+1: GDIM;
  GSIZET ibeg, iend; // beg, end indices for global array
  G_ .resize(nxy);
  G_ = NULLPTR;
  for ( auto j=0; j<nxy; j++ ) {
    G_ [j] = new StateComp(grid_->ndof());
  }

  Jac = &grid_->Jac();
  dXidX = &grid_->dXidX();


  // Cycle through all elements; fill metric elements
  for ( auto e=0; e<grid_->elems().size(); e++ ) {
    if ( (*gelems)[e]->elemtype() != GE_DEFORMED 
      && (*gelems)[e]->elemtype() != GE_2DEMBEDDED ) continue;

    ibeg = (*gelems)[e]->igbeg(); iend = (*gelems)[e]->igend();
    for ( auto j=0; j<GDIM; j++ ) {
      W[j]= (*gelems)[e]->gbasis(j)->getWeights();
      N[j]= (*gelems)[e]->size(j);
    }
    Jac->range(ibeg, iend);
    for ( auto j=0; j<nxy; j++ )
      for ( auto i=0; i<nxy; i++ ) (*dXidX)(i,j).range(ibeg, iend);

#if defined(_G_IS2D)

    for ( auto j=0; j<nxy; j++ ) { // G vector element 
      (*G_[j]).range(ibeg, iend); // restrict global vec to local range
      for ( auto m=0, n=0; m<N[1]; m++ ) {
        for ( auto l=0; l<N[0]; l++,n++ ) {
          (*G_[j])[n] = (*dXidX)(j,j)[n] 
                      * (*W[0])[l] * (*W[1])[m] * (*Jac)[n];
        }
      }
      (*G_[j]).range_reset(); // reset to global range
    }

#else

    for ( auto j=0; j<nxy; j++ ) { // G vector element 
      (*G_[j]).range(ibeg, iend); // restrict global vec to local range
        for ( auto p=0, n=0; p<N[2]; p++ ) {
          for ( auto m=0; m<N[1]; m++ ) {
            for ( auto l=0; l<N[0]; l++,n++ ) {
              (*G_[j])[n] = (*dXidX)(j,j)[n] 
                          * (*W[0])[l] * (*W[1])[m] * (*W[2])[p] * (*Jac)[n];
            }
          }
        }
      (*G_[j]).range_reset(); // reset global vec to global range
    }

#endif
  } // end, element list

  // Reset ranges to global scope:
  Jac->range_reset();
  for ( auto j=0; j<nxy; j++ )
    for ( auto i=0; i<nxy; i++ ) (*dXidX)(i,j).range_reset();


} // end of method def_init


//**********************************************************************************
//**********************************************************************************
// METHOD : reg_init
// DESC   : Compute metric components for regular elements. Mass *is not* included
//          here, so application must be done via mass operator.
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename TypePack>
void GpdV<TypePack>::reg_init()
{
  if ( grid_->itype(GE_REGULAR).size() <= 0 ) return; 


  GTVector<GSIZET>             N(GDIM);
  GTMatrix<StateComp>  *dXidX;    // element-based dXi/dX matrix
  StateComp            *Jac;      // element-based Jacobian
  GElemList                   *gelems = &grid_->elems();

  // Compute 'metric' components:
  // Gi = [dxi/dx, deta/dy, dzeta/dz]; 
  //
  // We have
  // 
  //    Jac = L1 L2 L3/8 and
  //  dXi_1/dx = 2 / L1
  //  dXi_2/dy = 2 / L2
  //  dXi_3/dz = 2 / L3
  //
  // We compute the metric terms as : Gii = (dXi_i/dX_i)^2
  // NOTE: unlike in the deformed case, we don't multiply
  //       by the weights or by Jacobian:
  //
  //  G1 = dXi^1/dX^1 = 2 / L1
  //  G2 = dXi^2/dX^2 = 2 / L2
  //  G3 = dXi^3/dX^3 = 2 / L3
  // These values should already be set in dXidX in grid

  Jac = &grid_->Jac();
  dXidX = &grid_->dXidX();

  GSIZET nxy = grid_->itype(GE_2DEMBEDDED).size() > 0 ? GDIM+1: GDIM;
  GSIZET ibeg, iend; // beg, end indices for global array
  G_ .resize(nxy);
  G_ = NULLPTR;
  for ( auto j=0; j<nxy; j++ ) {
    G_ [j] = new StateComp(grid_->ndof());
  }


  // Cycle through all elements; fill metric elements
  for ( auto e=0; e<grid_->elems().size(); e++ ) {
    if ( (*gelems)[e]->elemtype() != GE_REGULAR ) continue;

    for ( auto j=0; j<GDIM; j++ ) {
     *G_[j] = (*dXidX)(j,0);
      G_[j]->pointProd(*Jac);
    }
  } // end, element list


} // end of method reg_init


