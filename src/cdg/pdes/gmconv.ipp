//==================================================================================
// Module       : gmconv.ipp
// Date         : 6/11/20 (DLR)
// Description  : Object defining a moist convection solver:
//
//                PDEs:
//                  d_t rhoT + Div (rhoT v)   = -Ltot, total mass
//                  d_t q_v + v.Grad q_v      = q_v Ltot/rhoT + dot(s_v)/rhoT
//                  d_t q_h + v.Grad q_i      = q_h Ltot/rhoT - Div(rhoT q_h W_i)/rhoT
//                                            + dot(s_h)/rhoT
//
//                where 
//                  Ltot = Sum_h Div(rhoT q_h vec{W}_i), 
//                is the total mass loss due to hydrometeor fallout, and              
//                and q_h are the hydrometeor (liquid and ice) mass
//                fractions. The dry mass fraction is:
//                  q_d = 1 - Sum_h q_h.
//                Note:
//                  Sum_i dot(s_i) = 0, where sum is over all densities,
//                and dot(s_i) are the mass sources for vapor, and hydrometeors.
//                This solver can be built in 2D or 3D for box grids,
//                but is valid only for 3D spherical grids.
//
//                The State vector consists of the following:
//                  vx/sx   : x-velocity or momentum density 
//                  vy/sy   : y-velocity or momentum density 
//                  vz/sz   : z-velocity or momentum density 
//                  e       : internal or total energy density
//                  rho_tot : total density or dry density, if no moisture
//                  qvapor  : water vapor mass fraction
//                  q_liq_0 : liquid substance 0 mass fraction |
//                  q_liq_1 : liquid substance 1 mass fraction |  'liquid' mass sector
//                  q_liq_2 : liquid substance 2 mass fraction |
//                   ...
//                  q_ice_0 : 'ice' substance 0 mass fraction  |
//                  q_ice_1 : 'ice' substance 1 mass fraction  |  'ice' mass sector
//                  q_ice_2 : 'ice' substance 2 mass fraction  |
//                   ...
//                  w_liq_0 : liquid substance 2 term velocity |
//                  w_liq_1 : liquid substance 2 term velocity | 'liquid' term vel. sector
//                  w_liq_2 : liquid substance 2 term velocity |

//                  w_ice_0 : 'ice' substance 2 term velocity  |
//                  w_ice_1 : liquid substance 2 term velocity | 'ice' term vel. sector
//                  w_ice_2 : 'ice' substance 2 term velocity  |
//                   ...
//
//                The terminal velocities in this state are prescribed,
//                if used; all others are evolved. If hydrometeor fallout
//                is specified, then terminal velocities for all hydrometeors
//                must be provided. Each vector may be either of length 1, 
//                for which it's assumed to be a constant, or it may be 
//                space-dependent, for which it provides the terminal 
//                velocity in the preferred diection at all spatial locaions.
// 
// 
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : EquationBase.
//==================================================================================

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method  (1)
// DESC   : Instantiate with grid + state + tmp. 
//          grid      : grid object
//          traits    : GMConv:Traits struct
//          tmp       : Array of tmp vector pointers, pointing to vectors
//                      of same size as State. Must be MAX(2*DIM+2,iorder+1)
//                      vectors
// RETURNS: none
//**********************************************************************************
template<typename TypePack>
GMConv<TypePack>::GMConv(Grid &grid, GMConv<TypePack>::Traits &traits, State &tmp) :
EquationBase<TypePack>(),
bforced_                 (FALSE),
bsteptop_                (FALSE),
bvterm_                  (FALSE),
nevolve_                     (0),
nhydro_                      (0),
nmoist_                      (0),
gmass_                 (NULLPTR),
gimass_                (NULLPTR),
/*
gflux_                 (NULLPTR),
*/
ghelm_                 (NULLPTR),
gadvect_               (NULLPTR),
gexrk_                 (NULLPTR),
gpdv_                  (NULLPTR),
grid_                    (&grid),
comm_         (ggfx_->getComm()),
ggfx_         (&grid.get_ggfx()),
steptop_callback_      (NULLPTR)
{

  GINT      nexcl;
  GGridIcos *icos = dynamic_cast<GGridIcos*>(grid_);

  assert(tmp.size() >= req_tmp_size() && "Insufficient tmp space provided");
  assert(!(GDIM==2 && icos!=NULLPTR) && "Embedded 2D spherical grid not allowed");

  traits_.iforced.resize(traits.iforced.size());
  traits_ = traits;

  v_.resize(GDIM); v_ = NULLPTR;
  W_.resize(GDIM); W_ = NULLPTR;

  nexcl = 0;
  
  // Set space for state velocity, if needed:
  if ( traits_.usemomden ) {
    for ( auto j=0; j<GDIM; j++ ) {
      v_[j] = tmp[j];
      nexcl++;
    }
  }
  // Set space for terminal velocity, if needed:
  if ( icos != NULLPTR ) {
    for ( auto j=0; j<GDIM; j++ ) {
      W_[j] = tmp[j+v_.size()];
      nexcl++;
    }
  }

  // Set up tmp pool:
  utmp_.resize(tmp.size()-nexcl); 
  for ( auto j=nexcl; j<tmp.size(); j++ ) {
    utmp_[j] = tmp[j];
  }

  init();
  
} // end of constructor method (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
template<typename TypePack>
GMConv<TypePack>::~GMConv()
{
//if ( gflux_   != NULLPTR ) delete gflux_;
  if ( ghelm_   != NULLPTR ) delete ghelm_;
  if ( gadvect_ != NULLPTR ) delete gadvect_;
  if ( gpdv_    != NULLPTR ) delete gpdv_;
  if ( gexrk_   != NULLPTR ) delete gexrk_;

} // end, destructor


//**********************************************************************************
//**********************************************************************************
// METHOD : dt_impl
// DESC   : Compute time step, assuming a Courant number of 1:
//            dt = min_grid(dx/u)
// ARGS   : t : time
//          u : state
//          dt: timestep, returned
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::dt_impl(const Time &t, State &u, Time &dt)
{
   GString    serr = "GMConv<TypePack>::dt_impl: ";
   GINT       pmax;
   GSIZET     ibeg, iend;
   GFTYPE     dtmin, dt1, umax;
   GFTYPE     drmin  = grid_->minlength();
   GElemList *gelems = &grid_->elems();
   GTVector<GNBasis<GCTYPE,GFTYPE>*> *gbasis;

   // This is an estimate. We assume the timestep is
   // is governed by fast sonic waves with speed
   //  |v| + c,
   // where 
   //   c^2 = p/rho_dry; p ~ Rd e / Cv,
   // and e is int. energy density, and
   //   Cv = Cvd qd + Cvv qv + Sum_i(Cl_i ql_i) + Sum_j(Ci_j qi_j).
   // Then, dt is computed element-by-element from
   //   dt = dx_min/(v + c)_max
   // where min and max are computed over the element.
   
   dtmin = std::numeric_limits<GFTYPE>::max();

   // Compute Cv:
   utmp_[1] =           1.0;  
   utmp_[1] -= (*u[GDIM+1]);     // running total 1- Sum_k q_k
   utmp_[0] = (*u[GDIM+1])*CVV;  // Cv = Cvv * q_vapor
   for ( auto k=GDIM+3; k<traits_.nlsector+1; k++ ) { // vapor & liquids
     *utmp_[1] -= *u[k];         // subtract in ql_k
     *utmp_[0] += (*u[k]) * CVL; // add in Cvl * ql_k
   }
   for ( auto k=GDIM+3+traits_.nlsector; k<traits_.nisector; k++ ) { // ice
     *utmp_[1] -= *u[k];         // subtract in qi_k
     *utmp_[0] += (*u[k]) * CVI; // add in Cvi * qi_k
   }
   // After subtracting q_i, final result is qd=q_dry, so:
   *utmp_[0] += (*utmp_[1])*CVD; // Final Cv += Cvd * qd

   // Compute v^2:
   *utmp[1] = *u[0]; utmp[0]->rpow(2);
   for ( auto k=1; k<u.size(); k++ ) { // each advecting v 
     *utmp[2] = *u[k]; utmp[2]->rpow(2);
     *utmp[1] += (*utmp[2]);           // v^2 += v_k^2
   }
  
   // Compute sound speed: c^2 = p/rho_dry; p ~ Rd e / Cv,
   *utmp[2]  = (*u[GDIM]);   // utmp = e
   *utmp[2] /= *utmp[0];     // e / Cv
   *utmp[2] *= RD;           // c^2 = Rd e / Cv
   *utmp[2] += *utmp[1];     // v^2 + c^2
   GMTK::maxbyelem<GFTYPE>(*grid_, *utmp[2], maxbyelem_);
   
   // Note: maxbyelem_ is an array with the max of v^2 + c^2 
   //       on each element

   // Find estimate of smallest dt on this task:
   for ( auto e=1; e<dxmin_.size(); e++ ) { // check each element
     dt1 = dxmin_[e] / maxbyelem_[e]; // this dt^2
     dtmin = MIN(dtmin, sqrt(dt1)); 
   }

   // Find minimum dt over all tasks:
   GComm::Allreduce(&dtmin, &dt1, 1, T2GCDatatype<GFTYPE>() , GC_OP_MIN, comm_);

   // Limit any timestep-to-timestep increae to 10%:
   dt = MIN(dt1*traits_.courant, 1.1*dt);

} // end of method dt_impl


//**********************************************************************************
//**********************************************************************************
// METHOD : dudt_impl
// DESC   : Compute RHS for explicit schemes
// ARGS   : t   : time
//          u   : state. 
//          uf  : forcing tendency state vector
//          ub  : bdy state vector
//          dt  : time step
//          dudt: accelerations, returned
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::dudt_impl(const Time &t, const State &u, const State &uf,  const State &ub, const Time &dt, Derivative &dudt)
{
  assert(!traits_.bconserved &&
         "conservation not yet supported"); 

  GString    serr = "GMConv<TypePack>::dudt_impl: ";
  GINT       ibeg, nice, nliq ;
  StateComp *irhoT, *Ltot;
  StateComp *e, *p, *rhoT, *T; // energy den, pressure, temperature
  StateComp *Jac, *Mass;
  StateComp *tmp1, *tmp2;
  State      g(GDIM); 

  // NOTE:
  // Make sure that, in init(), Helmholtz op is using only
  // weak Laplacian (q * mass term isn't being used), or there will 
  // be problems. This is required for explicit schemes, for
  // which this method is called.

  assert( !traits_.bsonserved ); // don't allow conservative form yet

  // Assign qi, tvi, qice, qliq, tvice, tvliq:
  assign_helpers(u);

 *urhstmp_[0] = *u[GDIM+1]; urhstmp_[0]->rpow(-1.0); // 1/total mass

  // Set tmp pool for RHS computations:
  Ltot  = urhstmp_[urhstmp_.size()-1];
  irhoT = urhstmp_[urhstmp_.size()-2];
  tmp1  = urhstmp_[urhstmp_.size()-3];
  tmp2  = urhstmp_[urhstmp_.size()-4];

  // Compute velocity for timestep:
  compute_v(u, utmp_, v_); // stored in v_

  // Get 1/rhoT:
  rhoT  =  u[GDIM+1];
 *irhoT = *u[GDIM+1];
  irhoT->rpow(-1.0);
  

  // *************************************************************
  // Total density RHS:
  // *************************************************************

  compute_div(*u[GDIM+1], v_, urhstmp_, *dudt[GDIM+1]); 
  compute_falloutsrc(*u[GDIM+1], qi_, tvi_, -1, uoptmp_, *Ltot);
  GMTK::saxpby<Ftype>(*dudt[GDIM+1], 1.0, *Ltot, 1.0);   // += Ltot
  if ( uf[GDIM+1] != NULLPTR ) *dudt[GDIM+1] -= *uf[GDIM+1];//  += sdot(s_rhoT)
  
  // First, compute all operators as though they are on the LHS, then
  // change the sign and add Mass at the very end....

  // *************************************************************
  // Mass fraction equations (vapor + all hyrodmeteors) RHS:
  // We solve
  //  dq_i/dt + u.Grad q_i = -div(q_i rhoT W_i)/rhoT 
  //                       + q_i/rhoT Ltot + dot(s_i)/rhoT
  // where Ltot is total fallout source over liq + ice sectors:
  // *************************************************************
  ibeg   = GDIM + 2;
  for ( auto j=0; j<nmoist_; j++ ) {
    gadvect_->apply(*qi_[j], v_, uoptmp_, *dudt[ibeg+j]); // apply advection
    compute_vterm(*tvi_[j], W_);
    *tmp1 = (*qi_[j]) * (*rhoT);               // q_i rhoT
    compute_div(*tmp1, W_, urhstmp_, *tmp2);   // Div(q_i rhoT W)
    *tmp2 *= *irhoT;                           // Div(q_i rhoT W)/rhoT
    *dudt[ibeg+j] += *tmp2;                    // += Div(q_i rhoT W)/rhoT
    *tmp1  = (*Ltot) * (*irhoT); *tmp1 *= (*qi_[j]) // q_i/rhoT Ltot
    *dudt[ibeg+j] -= *tmp1;                    // += -q_i/rhoT Ltot
    if ( uf[ibeg+j] != NULLPTR ) {             // add in sdot(s_i)/rhoT
      *tmp1 = *uf[ibeg+1]; *tmp1 *= *irhoT;    // dot(s)/rhoT 
      *tmp1 *= *Jac ; *tmp1 *= *Mass;
      GMTK::saxpby<Ftype>(*dudt[ibeg+j], 1.0, *tmp1, -1.0); 
                                               // -= dot(s)/rhoT
    }
  }
  
  p = urhstmp_[urhstmp_.size()-1]; // holds pressure
  T = urhstmp_[urhstmp_.size()-2]; // holds temperature
  e = u[GDIM];                     // internal energy density

  // *************************************************************
  // Energy equationu RHS:
  // *************************************************************
  
  compute_ptemp(u, uoptmp_, *T, *p);               // find p, T

  GMTK::saxpby<Ftype>(*tmp1, *e, 1.0, *p, 1.0);    // h = p+e, enthalpy density
  compute_div(*tmp1, v_, uoptmp_, *dudt[GDIM]);    // Div (h v);

  if ( traits_.dofallout || traits_.dodry ) {
   *tmp1 = *rhoT; tmp1->pointProd(CPL, *T);        // tmp1 = CP_liq rhoT T
    compute_falloutsrc(*tmp1, qliq_, tvliq_, -1.0, uoptmp_, *Ltot);
                                                   // liquid fallout src
   *dudt[GDIM] += *Ltot;                           // += L_liq
   *tmp1 = *rhoT; tmp1->pointProd(CPI, *T);        // tmp1 = CP_ice rhoT T
    compute_falloutsrc(*tmp1, qice_, tvice_, -1.0, uoptmp_, *Ltot); 
                                                   // ice fallout src
   *dudt[GDIM] += *Ltot;                           // += L_ice
  }

  gadvect_->apply(*p, v_, uoptmp_, *tmp1);         // v.Grad p
 *dudt[GDIM] -= *tmp1;                             // -= v . Grad p

  if ( traits_.dograv ) {
    compute_pe(*rhoT, qi_, tvi_, uoptmp_, *tmp1);
   *dudt[GDIM] += *tmp1;                           // += Sum_i rhoT q_i g.W_i
  }

  compute_fv(fv_, v_, *tmp2, *tmp1);
 *dudt[GDIM] += *tmp1;                             // += f_kinetic . v

  if ( uf[GDIM] != NULLPTR ) {                    
    *tmp1 = *uf[GDIM]; *tmp1 *= *Jac ; *tmp1 *= *Mass;
    GMTK::saxpby<Ftype>(*dudt[GDIM], 1.0, *tmp1, -1.0); 
                                                   // -= q_heat
  }

  // *************************************************************
  // Momentum equations RHS:
  // *************************************************************
  Jac = &grid_->Jac();
  Mass = grid_->massop().data();
  if ( traits_.dograv) {
    *tmp1 = -GG; // set constant grav accel
     for ( auto j=0; j<GDIM; j++ ) g[j] = utmp[j];
     compute_vterm(*tmp1, g); // compute grav vector components
  }

  for ( auto j=0; k<GDIM; j++ ) {

    gadvect_->apply(*s_[j], v_, uoptmp_, *dudt[j]);   // v.Grad s_j
    if ( traits_.dofallout || traits_.dodry ) {
      compute_falloutsrc(*u[j], qli_, tvli_,-1.0, uoptmp_, *Ltot);
                                                      // hydrometeor fallout src
     *dudt[GDIM] += *Ltot;                            // += L_tot
    }
    grid_->wderiv(*u[j], j+1, TRUE, uoptmp_, *tmp1);  // Grad p
   *dudt[GDIM] += *tmp1;                              // += Grad p
    ghelm_->opVec_prod(*u[j], uoptmp_, *tmp1);        // nu Laplacian s_j
   *dudt[GDIM] -= *tmp1;                              // -= nu Laplacian s_j
    if ( traits_.docoriolis ) {
      GMTK::cross_prod(g, s_, j+1, *tmp1);
     *tmp1 *= *Jac; *tmp1 *= *Mass;             
      GMTK::saxpby<Ftype>(*dudt[j], 1.0, *tmp1, 2.0); // += 2 Omega X (thooT v) M J
    }
    if ( traits_.dograv && g[j] != NULLPTR ) {
      g[j]->pointProd(*rhotT, *tmp1);
     *tmp1 *= *Jac; *tmp1 *= *Mass;             
     *dudt[GDIM] -= *tmp1;                            // -= rhoT g M J
    }
    if ( traits_.forced && uf[j] != NULLPTR ) {                    
      *tmp1 = *uf[j]; *tmp1 *= *Jac ; *tmp1 *= *Mass;
      GMTK::saxpby<Ftype>(*dudt[j], 1.0, *tmp1, -1.0); 
                                                      // -= f_v
    }

  } // end, momentum loop

  // Multiply RHS by -M^-1 to (1) place all terms on the 'RHS',
  // and (2) to account for factor of M on du/dt term:
  for ( auto j=0; j<nevolve_; j++ ) {
    dudt[j]->pointProd(-1.0, *gimass_->data());// dudt -> -M^-1 dudt
  }
  
} // end of method dudt_impl


//**********************************************************************************
//**********************************************************************************
// METHOD : step_impl (1)
// DESC   : Step implementation method entry point
// ARGS   : t   : time
//          uin : input state, modified on output with update
//          uf  : force-tendency vector
//          ub  : bdy vector
//          dt  : time step
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::step_impl(const Time &t, State &uin, State &uf, State &ub, const Time &dt)
{

  GBOOL bret;

  // If there's a top-of-the-timestep callback, 
  // call it here:
  if ( bsteptop_ ) {
    steptop_callback_(t, uin, dt);
  }


  // Set evolved state vars from input state.
  // These are not deep copies:
  for ( auto j=0; j<uin.size(); j++ ) uevolve_ [j] = uin[j];

  switch ( traits_.isteptype ) {
    case GSTEPPER_EXRK:
      for ( auto j=0; j<uold_.size(); j++ ) *uold_[j] = *uevolve_[j];
      step_exrk(t, uold_, uf, ub, dt, uevolve_);
      break;
    case GSTEPPER_BDFAB:
    case GSTEPPER_BDFEXT:
      step_multistep(t, uevolve_, uf, ub, dt);
      break;
  }

  // Check solution for NaN and Inf:
  bret = TRUE;
  for ( auto j=0; j<uevolve_.size(); j++ ) {
     bret = bret && uevolve_ [j]->isfinite();
  }
  assert(bret && "Solution not finite");

} // end of method step_impl (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : step_impl (2)
// DESC   : Step implementation method entry point
// ARGS   : t   : time
//          uin : input state, modified on output with update
//          ub  : bdy vector
//          dt  : time step
//          uout: output state
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::step_impl(const Time &t, const State &uin, State &uf,  State &ub, const Time &dt, State &uout)
{
  assert(FALSE && "step_impl(2) not available");

} // end of method step_impl (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : step_multistep
// DESC   : Carries out multistep update. The time derivative and 
//          advection terms are handlex using a multistep expansion. The
//          advection term is 'extrapolated' to the new time level
//          using known state data so as to obviate the need for a fully 
//          implicit treatment. The dissipation term is handled implicitly.
// ARGS   : t   : time
//          u   : state
//          uf  : force tendency vector
//          ub  : bdy vector
//          dt  : time step
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::step_multistep(const Time &t, State &uin, State &uf, State &ub, const Time &dt)
{
  assert(FALSE && "Multistep methods not yet available");
  
} // end of method step_multistep


//**********************************************************************************
//**********************************************************************************
// METHOD : step_exrk
// DESC   : Take a step using Explicit RK method
// ARGS   : t   : time
//          uin : input state; must not be modified
//          uf  : force tendency vector
//          ub  : bdy vector
//          dt  : time step
//          uout: output/updated state
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::step_exrk(const Time &t, State &uin, State &uf, State &ub, const Time &dt, State &uout)
{
  assert(gexrk_ != NULLPTR && "GExRK operator not instantiated");


  // GExRK stepper steps entire state over one dt:
  gexrk_->step(t, uin, uf, ub, dt, urktmp_, uout);

//GMTK::constrain2sphere(*grid_, uout);

} // end of method step_exrk


//**********************************************************************************
//**********************************************************************************
// METHOD : init
// DESC   : Initialize equation object
// ARGS   : traits: GMConv::Traits variable
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::init()
{
  GString serr = "GMConv<TypePack>::init: ";

  GBOOL      bmultilevel = FALSE;
  GSIZET     n;
  GSIZET     nc = grid_->gtype() == GE_2DEMBEDDED ? 3 : GDIM;
  GINT       nop;
  CompDesc *icomptype = &this->stateinfo().icomptype;

  // Check if specified stepper type is valid:
  GBOOL  bfound;
  GSIZET itype;
  valid_types_.push_back("GSTEPPER_EXRK");
/*
  valid_types_.push_back("GSTEPPER_BDFAB");
  valid_types_.push_back("GSTEPPER_BDFEXT");
*/
  bfound = valid_types_.contains(traits_.ssteptype, itype);
  assert( bfound && "Invalid stepping method specified");
  traits_.isteptype = static_cast<GStepperType>(itype);

  if ( traits_.dodry ) traits_.dofallout = FALSE;

  // Set dissipation from traits. Note that
  // this is set to be constant, based on configuration,
  // even though the solver can accommodate spatially
  // variable dissipation:
  nu_.resize(1);
  nu_ = traits_.nu;

  // Find no. state and solve members, and component types:
  for( auto j=0; j<GDIM; j++ ) icomptype->push_back(GSC_KINETIC); 
  icomptype->push_back(GSC_DENSITYT); 
  icomptype->push_back(GSC_ENERGY); 
  if ( traits_.dodry )  {
    traits_.nsolve = GDIM + 2;
    traits_.nstate = traits_.nsolve;
  }
  else {
    n = traits_.nlsector + traits_.nisector;
    traits_.nsolve = GDIM + 2 + n;
    traits_.nstate = traits_.nsolve;
    for( auto j=0; j<traits_.nlsector; j++ ) icomptype->push_back(GSC_MASSFRAC); 
    for( auto j=0; j<traits_.nisector; j++ ) icomptype->push_back(GSC_MASSFRAC); 
    if ( traits_.dofallout ) { // include terminal velocities:
      traits_.nstate += 
                     traits_.nlsector + traits_.nisector;
    for( auto j=0; j<n; j++ ) icomptype->push_back(GSC_PRESCRIBED); 
    }
  }

  // Find multistep/multistage time stepping coefficients:
  GMultilevel_coeffs_base<GFTYPE> *tcoeff_obj=NULLPTR; // time deriv coeffs
  GMultilevel_coeffs_base<GFTYPE> *acoeff_obj=NULLPTR; // adv op. coeffs


  std::function<void(const Time &t,                    // RHS callback function
                     const State  &uin,
                     const State  &uf,
                     const State  &ub,
                     const Time   &dt,
                     State &dudt)> rhs
                  = [this](const Time &t,           
                     const State  &uin, 
                     const State  &uf, 
                     const State  &ub, 
                     const Time   &dt,
                     State &dudt){dudt_impl(t, uin, uf, ub, dt, dudt);}; 

  std::function<void(const Time &t,                    // Bdy apply callback function
                     State  &uin, 
                     State &ub)> applybc 
                  = [this](const Time &t,              
                     State  &uin, 
                     State &ub){apply_bc_impl(t, uin, ub);}; 

  // Configure time stepping:
  switch ( traits_.isteptype ) {
    case GSTEPPER_EXRK:
      gexrk_ = new GExRKStepper<GFTYPE>(*grid_, traits_.itorder);
      gexrk_->setRHSfunction(rhs);
      gexrk_->set_apply_bdy_callback(applybc);
      gexrk_->set_ggfx(ggfx_);
      // Set 'helper' tmp arrays from main one, utmp_, so that
      // we're sure there's no overlap:
      uold_   .resize(traits_.nsolve); // solution at time level n
      urktmp_ .resize(traits_.nsolve*(traits_.itorder+1)+1); // RK stepping work space
      urhstmp_.resize(1); // work space for RHS
      nop = utmp_.size()-uold_.size()-urktmp_.size()-urhstmp_.size();
      assert(nop > 0 && "Invalid operation tmp array specification");
      uoptmp_ .resize(nop); // RHS operator work space
      // Make sure there is no overlap between tmp arrays:
      n = 0;
      for ( GSIZET j=0; j<traits_.nsolve ; j++, n++ ) uold_   [j] = utmp_[n];
      for ( GSIZET j=0; j<urktmp_ .size(); j++, n++ ) urktmp_ [j] = utmp_[n];
      for ( GSIZET j=0; j<urhstmp_.size(); j++, n++ ) urhstmp_[j] = utmp_[n];
      for ( GSIZET j=0; j<uoptmp_ .size(); j++, n++ ) uoptmp_ [j] = utmp_[n];
      break;
/*
    case GSTEPPER_BDFAB:
      dthist_.resize(MAX(traits_.itorder,traits_.inorder));
      tcoeff_obj = new G_BDF<GFTYPE>(traits_.itorder, dthist_);
      acoeff_obj = new G_AB<GFTYPE> (traits_.inorder, dthist_);
      tcoeffs_.resize(tcoeff_obj->getCoeffs().size());
      acoeffs_.resize(acoeff_obj->getCoeffs().size());
      tcoeffs_ = tcoeff_obj->getCoeffs(); 
      acoeffs_ = acoeff_obj->getCoeffs();
      uold_   .resize(traits_.nsolve); // solution at time level n
      for ( GSIZET j=0; j<traits_.nsolve; j++ ) uold_[j] = utmp_[j];
      bmultilevel = TRUE;
      break;
    case GSTEPPER_BDFEXT:
      dthist_.resize(MAX(traits_.itorder,traits_.inorder));
      tcoeff_obj = new G_BDF<GFTYPE>(traits_.itorder, dthist_);
      acoeff_obj = new G_EXT<GFTYPE>(traits_.inorder, dthist_);
      tcoeffs_.resize(tcoeff_obj->getCoeffs().size());
      acoeffs_.resize(acoeff_obj->getCoeffs().size());
      tcoeffs_ = tcoeff_obj->getCoeffs(); 
      acoeffs_ = acoeff_obj->getCoeffs();
      urhstmp_.resize(utmp_.size()-urktmp_.size());
      for ( GSIZET j=0; j<utmp_.size(); j++ ) urhstmp_[j] = utmp_[j];
      bmultilevel = TRUE;
      break;
*/
    default:
      assert(FALSE && "Invalid stepper type");
  }
  if ( tcoeff_obj != NULLPTR ) delete tcoeff_obj;
  if ( acoeff_obj != NULLPTR ) delete acoeff_obj;
  
  // Instantiate spatial discretization operators:
  gmass_   = &grid_->massop();
  ghelm_   = new GHelmholtz(*grid_);

  ghelm_->set_Lap_scalar(nu_);

  if ( traits_.isteptype ==  GSTEPPER_EXRK ) {
    gimass_ = &grid_->imassop();
  }

  // If doing semi-implicit time stepping; handle viscous term 
  // (linear) inplicitly, which implies using full Helmholtz operator:
  if ( traits_.isteptype == GSTEPPER_BDFAB || traits_.isteptype == GSTEPPER_BDFEXT ) {
    assert(FALSE && "Implicit time stepping not yet supported");
  }

  if ( traits_.bconserved ) {
    assert(FALSE && "Conservation not yet supported");
    gpdv_  = new GpdV<TypePack>(*grid_);
//  gflux_ = new GFlux(*grid_);
    assert( (gmass_   != NULLPTR
          && ghelm_   != NULLPTR
          && gpdv_    != NULLPTR) && "1 or more operators undefined");
  }
  else {
    gadvect_ = new GAdvect(*grid_);
    gpdv_    = new GpdV<TypePack>(*grid_);
    assert( (gmass_   != NULLPTR
          && ghelm_   != NULLPTR
          && gpdv_    != NULLPTR
          && gadvect_ != NULLPTR) && "1 or more operators undefined");
  }

  // If doing a multi-step method, instantiate (deep) space for 
  // required time levels for state:
  if ( bmultilevel ) {
    ukeep_ .resize(traits_.itorder);
    for ( auto i=0; i<traits_.itorder-1; i++ ) { // for each time level
      ukeep_[i].resize(traits_.nsolve);
      for ( auto j=0; j<traits_.nsolve; j++ ) ukeep_[i][j] = new GTVector<GFTYPE>(grid_->ndof());
    }
  }

  // Find minimum element edge/face lengths for 
  // timestep computation:
  maxbyelem_.resize(grid_->nelems());
  grid_->minlength(&dxmin_);

  // Set size of mass frac and term vel vectors,
  // misc. helper arrays:
  nhydro_ = traits_.dodry ? 0 : traits_.nlsector + traits_.nisector + 1;
  nmoist_ = traits_.dodry ? 0 : nhydro_ + 1;
  nevolve_ = GDIM + 2 + nmoist_;
  this->stateinfo.nevolve = traits_.solve;
  this->stateinfo.presc   = traits_.nstate - traits_.solve;
  qi_   .resize(nmoist_);
  tvi_  .resize(nmoist_);
  qliq_ .resize(traits_.nlsector);
  tvliq_.resize(traits_.nlsector);
  qice_ .resize(traits_.nisector);
  tvice_.resize(traits_.nisector);
  fk_   .resize(GDIM); 
  s_   . resize(GDIM); 

  qi_   = NULLPTR;
  tvi_  = NULLPTR;
  qliq_ = NULLPTR;
  tvliq_= NULLPTR;
  qice_ = NULLPTR;
  tvice_= NULLPTR;
  fk_   = NULLPTR;
  s_    = NULLPTR;

} // end of method init


//**********************************************************************************
//**********************************************************************************
// METHOD : cycle_keep
// DESC   : Cycle the mult-level states making sure the most
//          recent is at index 0, the next most recent, at index 1, etc...
// ARGS   : u     : State variable providing most recent state
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::cycle_keep(State &u)
{

  // Make sure following index map contains the 
  // correct time level information:
  //   ukeep[0] <--> time level n (most recent)
  //   ukeep[1] <--> time level n-1
  //   ukeep[2] <--> time level n-2 ...
  ukeep_ .resize(traits_.itorder);
  for ( auto i=traits_.itorder-1; i>=1; i-- ) ukeep_[i] = ukeep_[i+1];
  ukeep_[0] = u;

} // end of method cycle_keep


#if 0
//**********************************************************************************
//**********************************************************************************
// METHOD : set_nu
// DESC   : Set viscosity quantity. This may be a field, or effectively a
//          scalar (where only element 0 contains valid data). If not set, 
//          Helmholtz op creates a default of nu = 1.
//          nu : viscosity vector (may be of length 1). Will be managed
//               by caller; only a pointer is used by internal operators.
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::set_nu(GTVector<GFTYPE> &nu)
{
  assert(ghelm_ != NULLPTR && "Init must be called first");
  nu_ = nu; // Not sure this class actually needs this. May be removed later
  ghelm_->set_Lap_scalar(nu_);

} // end of method set_nu

#endif


//**********************************************************************************
//**********************************************************************************
// METHOD : apply_bc_impl
// DESC   : Apply global domain boundary conditions, ub
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::apply_bc_impl(const Time &t, State &u, State &ub)
{
  Time ttime = t;

  BdyUpdateList *updatelist = &grid_->bdy_update_list();;


  // Update bdy values if required to:
  for ( auto k=0; k<updatelist->size(); k++ ) { // foreach grid bdy
    for ( auto j=0; j<(*updatelist)[j].size(); j++ ) { // each update method
      (*updatelist)[k][j]->update(*grid_, this->stateinfo_, ttime, utmp_, u, ub);
    }
  }

} // end of method apply_bc_impl


//**********************************************************************************
//**********************************************************************************
// METHOD : req_tmp_size
// DESC   : Find required tmp size on GLL grid
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
GINT GMConv<TypePack>::req_tmp_size()
{
  GINT isize = 0;
 
  isize  = 2*GDIM + 3;

  if ( traits_.isteptype == GSTEPPER_EXRK ) {
    isize += traits_.nstate * traits_.itorder; 
  }

  return isize;
  
} // end of method req_tmp_size


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_cv
// DESC   : Compute total specific heat at const vol, Cv
//             Cv = Cvd qd + Cvv qv + Sum_i(Cl_i ql_i) + Sum_j(Ci_j qi_j).
//          where ql are the liquid mass fractions, and qi are the ice
//          mass fractions. 
// ARGS   : u    : state
//          utmp : tmp vectors; 1 required; only first one used
//          cv   : cv field
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::compute_cv(State &u, State &utmp, StateComp &cv)
{
   GString    serr = "GMConv<TypePack>::compute_cv: ";
   GINT       ibeg;

   // Compute Cv:
   if ( traits_.dodry ) { // if dry dynamics only
     utmp[0] = CVD;  
     return;
   }

   utmp[0] =           1.0;  
   utmp[0] -= (*u[GDIM+1]);      // running total 1- Sum_k q_k
   cv       = (*u[GDIM+1])*CVV;  // Cv = Cvv * q_vapor
   ibeg = GDIM + 3;
   for ( auto k=ibeg; k<ibeg+traits_.nlsector+1; k++ ) { // liquids
     *utmp[0] -= *u[k];         // subtract in ql_k
     *cv      += (*u[k]) * CVL; // add in Cvl * ql_k
   }
   ibeg = GDIM + 3 + traits_.nlsector;;
   for ( auto k=ibeg; k<ibeg+traits_.nisector; k++ ) { // ices
     *utmp[0] -= *u[k];         // subtract in qi_k
     cv       += (*u[k]) * CVI; // add in Cvi * qi_k
   }
   // After subtracting q_i, final result is qd=q_dry, so:
   cv  += (*utmp_[0])*CVD; // Final Cv += Cvd * qd

} // end of method compute_cv


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_qd
// DESC   : Compute dry mass fraction from other mass fractions:
//             Cv = 1 - Sum_i q_i
// ARGS   : u    : state
//          utmp : tmp vectors; at least 1 required; only first one used
//          qd   : dry mass fraction field
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::compute_qd(State &u, State &utmp, StateComp &qd)
{
   GString    serr = "GMConv<TypePack>::compute_qd: ";
   GINT       ibeg;

   assert(utmp.size() >= 1);

   // Compute qd:
   if ( traits_.dodry ) { // if dry dynamics only
     qd = 1.0;  
     return;
   }

   qd  =           1.0;  
   qd -= (*u[GDIM+1]);      // running total 1- Sum_k q_k
   ibeg = GDIM + 3;
   for ( auto k=ibeg; k<ibeg+traits_.nlsector+1; k++ ) { // liquids
     qd -= *u[k];         // subtract in ql_k
   }
   ibeg = GDIM + 3 + traits_.nlsector;;
   for ( auto k=ibeg; k<ibeg+traits_.nisector; k++ ) { // ices
     qd  -= *u[k];         // subtract in qi_k
   }

} // end of method compute_qd


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_temp
// DESC   : Compute temperature from state
//             T = eps / Cv = e_s /( rho' * Cv ),
//          with e_s the sensible internal energy density, 
//          rho' = total density fluctuations,
//             Cv = Cvd qd + Cvv qv + Sum_i(Cl_i ql_i) + Sum_j(Ci_j qi_j).
// ARGS   : u    : state
//          utmp : tmp vectors; 2 required; only first two used
//          temp : temperature field
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::compute_temp(State &u, State &utmp, StateComp &temp)
{
   GString    serr = "GMConv<TypePack>::compute_temp: ";
   StateComp *cv, *d, *e; 

   assert(utmp.size() >= 2);

   // Set int energy and density:
   e  = u[GDIM];   // sensible internal energy density
   d  = u[GDIM+1]; // total density fluctuation
   cv = utmp[1];   // Cv

   // Get Cv:
   compute_cv(u, utmp, *cv); // utmp[0] only used in call

   // Compute temperature:
   for ( auto j=0; j<es->size(); j++ ) {
     temp[j] = (*e)[j] / ( (*d)[j] * (*cv)[j] );
   }

} // end of method compute_temp


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_p 
// DESC   : Compute total pressure fluctuations from state
//              p' = rho' ( qd Rd + qv Rv ) T,
//          with total density fluctuation, rho', qd, qv the
//          dry and vapor mass fractions, Rd, Rv, the dry and vapor
//          gas constants, and T the temperature.
// ARGS   : u    : state
//          utmp : tmp vectors; 3 required; only first 3 used
//          p    : pressure fluctuation field
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::compute_p(State &u, State &utmp, StateComp &p)
{
   GString    serr = "GMConv<TypePack>::compute_p: ";
   StateComp *qd, *qv, *t; 


   assert(utmp.size() >= 3);

   // Set int energy and density:
   es = u[GDIM];   // sensible internal energy density
   d  = u[GDIM+1]; // total density fluctuation

   t = utmp[2];    // temp
   compute_temp(u, utmp, *t);  // first 2 utmp arrays used
  
   if ( traits_.dodry ) { // if dry dynamics only
     // p' = rho'  Rd T:
     for ( auto j=0; j<p.size(); j++ ) {
       p[j] = (*d)[j] * RD  * (*t)[j];
     }
     return;
   }

   // Set vapor mass fraction:
   qv = u[GDIM+2]; // qv, from state

   // Get dry mass fraction:
   qd = utmp[1];
   compute_qd(u, utmp, *qd); // first utmp array used 

   // p' = rho' ( qd Rd + qv Rv ) T:
   for ( auto j=0; j<p.size(); j++ ) {
     p[j] = (*d)[j] * ( (*qd)[j]*RD + (*qv)[j]*RV ) * (*t)[j];
   }

} // end of method compute_p 


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_ptemp 
// DESC   : Compute total pressure & temperature fluctuations from state
//              p' = rho' ( qd Rd + qv Rv ) T,
//          with total density fluctuation, rho', qd, qv the
//          dry and vapor mass fractions, Rd, Rv, the dry and vapor
//          gas constants, and T the temperature:
//             T = eps / Cv = e_s /( rho' * Cv ),
//          where e_s is the sensible internal energy density
// ARGS   : u    : state
//          utmp : tmp vectors; 2 required; only first 2 used
//          temp : array to hold temperature
//          p    : pressure fluctuation field
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::compute_ptemp(State &u, State &utmp, StateComp &temp, StateComp &p)
{
   GString    serr = "GMConv<TypePack>::compute_ptemp: ";
   StateComp *qd, *qv, *t; 


   assert(utmp.size() >= 3);

   // Set int energy and density:
   es = u[GDIM];   // sensible internal energy density
   d  = u[GDIM+1]; // total density fluctuation

   t = &temp;    // temp
   compute_temp(u, utmp, *t);  // first 2 utmp arrays used
  
   if ( traits_.dodry ) { // if dry dynamics only
     // p' = rho'  Rd T:
     for ( auto j=0; j<p.size(); j++ ) {
       p[j] = (*d)[j] * RD  * (*t)[j];
     }
     return;
   }

   // Set vapor mass fraction:
   qv = u[GDIM+2]; // qv, from state

   // Get dry mass fraction:
   qd = utmp[1];
   compute_qd(u, utmp, *qd); // first utmp array used 

   // p' = rho' ( qd Rd + qv Rv ) T:
   for ( auto j=0; j<p.size(); j++ ) {
     p[j] = (*d)[j] * ( (*qd)[j]*RD + (*qv)[j]*RV ) * (*t)[j];
   }

} // end of method compute_ptemp


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_div
// DESC   : Compute flux divergence  of quantity q:
//             Div ( q v)
//          with the intent of handling this conservatively or not.  
// ARGS   : q    : quantiy whose flux we compute divergence of
//          v    : vector of velocity components
//          utmp : tmp vectors; GDIM+1 required; only first GDIM+1 used
//          div  : divergence result
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::compute_div(StateComp &q, State &v, State &utmp, StateComp &div)
{
   GString    serr = "GMConv<TypePack>::compute_div: ";
   State     *tmp(GDIM);
   StateComp *qd, *qv, *t; 

   assert(utmp.size() >= GDIM+1);

   if ( traits_.bconserved ) {
     assert(FALSE); // conserved form not available yet
   }
   else {
     //   Div (q v) = q Div v + v.Grad q 
     assert(gadvect_ != NULLPTR && gpdv_ != NULLPTR);    
     for ( auto j=0; j<GDIM; j++ ) tmp[j] = utmp[j];
     gadvect->apply(q, v, tmp, div); 
     gpdv   ->apply(q, v, tmp, *utmp[GDIM]); 
     div += *utmp[GDIM];
   }

} // end of method compute_div


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_v
// DESC   : Compute velocity from momentum density in state vector.
//             v_i = s_i/rhot,
//          where v_i is the member data array.
// ARGS   : u    : state
//          utmp : tmp vectors; first array used. If traits.usemomden==FALSE,
//                 then no tmp arrays are required
//          v    : velocity state; components may change on exit
// RETURNS: none. Member data, v_, is set here
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::compute_v(State &u, State &utmp, State &v)
{
   GString    serr = "GMConv<TypePack>::compute_v: ";
   State     *tmp(GDIM);
   StateComp *irhot, *rhot

   assert(utmp.size() >= 1);

   if ( !traits_.usemomden ) {
     // State uses velocity form already so 
     // do pointer assignment:
     for ( auto j=0; j<GDIM; j++ ) { // v_i = u_i
       v[j] = u[j];
     }
     return;
   }

   // Compute inverse mass:
   rhot  = u[GDIM+1]; 
   irhot = utmp[0];
   for ( auto j=0; j<rhot->size(); j++ ) (*irhot)[j] = 1.0/(*rhot)[j];

   // Find velocity from momentum density:
   
   for ( auto j=0; j<GDIM; j++ ) {
      v[j]  = v_[j];    // set to allocated member data
     *v[j]  = *u[j];    // deep copy
     *v[j] *= (*irhot); // divide by density
   }

} // end of method compute_v


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_vterm
// DESC   : Transform 'terminal' (preferred) direction component specified into
//          Cartesian compoments, W.
// ARGS   : tvi   : terminal velocity or gravity magnitude in preferred
//                  direction
//          W     : Cartesian vector field. Components may be set to NULL if not
//                  needed.
// RETURNS: none. 
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::compute_vterm(StateComp &tvi, State &W)
{
   GString    serr = "GMConv<TypePack>::compute_vterm: ";
   GINT       ibeg, nhydro;
   Ftype      r, x, y, z;
   Ftype      lat, long;
   GTVector<GTVector<GFTYPE>> 
             *xnodes = &grid.xNodes();

   GGridIcos *icos = dynamic_cast<GGridIcos*>(grid_);
   GGridIcos *box  = dynamic_cast <GGridBox*>(grid_);

   assert(W.size() == GDIM);

   // Specify components for Cartesian grid:
   if ( box != NULLPTR ) {
     // In Cartesian coords, select the 'z' direction
     // as preferred 'fallout' direction. In 2d, this
     // will be the 2-coord; in 3d, the 3-coord:
     W = NULLPTR;
     W[GDIM-1] = tvi; 
     return;
   }

   // Specify components for spherical grid. Here,
   // the preferred direction is along radial direction; 
   // we assume that the prescribed term velocity
   // is the radial component only, which must be
   // transformed into Cartesian components:

   for ( auto i=0; i<(*xnodes)[0].size(); i++ ) {
     x = (*xnodes)[0][i]; y = (*xnodes)[1][i]; z = (*xnodes)[2][i];
     r     = sqrt(x*x + y*y + z*z);
     lat   = asin(z/r);
     long  = atan2(y,x);

     (*W[0])[i] = cos(lat)*cos(long);
     (*W[1])[i] = cos(lat)*sin(long);
     (*W[2])[i] = sin(lat);
   }

   for ( auto j=0; j<W_.size(); j++ ) {
     *W[j] *= (*tvi);
   }


} // end of method compute_vterm


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_falloutsrc
// DESC   : Compute total effect due to fallout of density quantity, g. 
//            r = Sum_i Div [q_i g vec{W_i}] 
//          where vec{W_i} is the fallout terminal velocity for precipitating
//          species i (from either liquid or ice sectors of state). 
// ARGS   : 
//          g    : quantiy to fall out ('flux-out'). May be: total density,
//                 momentum density component, internal energy density
//          qi   : hydrometeor mass fractions
//          tvi  : terminal velocity vector for each qi hydrometeor, from which
//                 vec{W} is computed
//          jexcl: which hydrometeor index to exclude from sum. If jexcl<0,
//                 exclude none. This index must be the 0-starting index of
//                 the hydrometeor in the qi array.
//          utmp : tmp vectors; GDIM+3 required
//          r    : fallout src field
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::compute_falloutsrc(StateComp &g, State &qi, State &tvi, GINT jexcl, State &utmp, StateComp &r)
{
   GString     serr = "GMConv<TypePack>::compute_falloutsrc: ";
   GINT        nhydro; // no. hydrometeors
   StateComp  *div, *qg; 

   assert(tvi.size() == qi.size());
   assert(utmp.size() >= 2*GDIM+1);
   assert(jexcl < 0 || (jexcl >=0 && jexcl < qi.size()));

   // Compute:
   //    r = -Sum_i Div (rho_t q_i vec{W}_i )

   r = 0.0;
   if ( !traits_.dofallout || traits_.dodry ) return;

   nhydro = qi.size(); // traits_.nlsector + traits_.nisector;

   qg    = utmp[GDIM+1];    // temp
   div   = utmp[GDIM+2];    // temp
   for ( auto j=0; j<nhydro; j++ ) {
     if ( j == jexcl ) continue;  

     *qg = g; (*qg) *= (*qi[j]); // compute g q_i

     // Convert terminal velocities to required 
     // (Cartesian) components. vterm may point
     // to tvi, or to v_ member data, or may be NULL:
     compute_vterm(*tvi[j], W_);

     // Compute i_th contribution to source term:
     compute_div(*qg, W_, utmp, *div);
     r += *div;
   }

} // end of method compute_falloutsrc


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_fv
// DESC   : Compute dot prod of input vector fields,
//          with the understanding that either vector 
//          could have NULL or constant components:
//                r = f . v,
//          where return vector, r, is non-NULL and of
//          full length. A constant vector has a size of 1
// ARGS   : 
//          f  : forcing vector; may be NULL or constant
//          v  : velocity components; may be NULL or constant
//          tmp: tmp vector, of full size
//          r  : scalar field containing injection rate, of full size
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::compute_fv(State &f, State &v, StateComp &tmp, StateComp &r)
{
   GString     serr = "GMConv<TypePack>::compute_fv: ";

   assert(f.size() == v.size());


   r = 0.0;
   for ( auto j=0; j<f.size(); j++ ) {

     if ( f[j] == NULLPTR || v[j] == NULLPTR ) continue;
     if      ( f[j]->size() == 1 && v[j]->size() == 1 ) { // f , v constant
       tmp =  (*f[j])[0] * (*v[j])[0];
     }
     else if ( f[j]->size() >  1 && v[j]->size() == 1 ) { // v constant
       GMTK::saxpby<Ftype>(tmp, 1.0, *f[j], (*v[j])[0]);  
       
     }
     else if ( f[j]->size() == 1 && v[j]->size() >  1 ) { // f constant
       GMTK::saxpby<Ftype>(tmp, 1.0, *v[j], (*f[j])[0]);  
     }
     else {                                               // f, v of full length
       f[j]->pointProd(*v[j], tmp);
     }
     r += tmp;

   }

} // end of method compute_fv


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_pe
// DESC   : Compute energy potential energy of hydrometeors:
//       
//              r = Sum_i rhoT q_i vec{g} . vec{W}_i,
//          where i is over all hydrometeor sectors, q_i & W_i
//          are mass fraction and terminal velocity of ith hydrometeor. 
// ARGS   : 
//          rhoT : total density
//          qi   : hydrometeor mass fractions
//          tvi  : terminal velocity vector for each qi hydrometeor, from which
//                 vec{W} is computed
//          utmp : tmp vectors; GDIM+2 required
//          r    : fallout potential energy
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::compute_pe(StateComp &rhoT, State &qi, State &tvi, State &utmp, StateComp &r)
{
   GString     serr = "GMConv<TypePack>::compute_pe: ";
   GINT        nhydro; // no. hydrometeors
   StateComp   *tmp1, *tmp2;
   State       g(GDIM); 

   assert(tvi.size() == qi.size());
   assert(utmp.size() >= GDIM+2);
   assert(jexcl < 0 || (jexcl >=0 && jexcl < qi.size()));

   // Compute:
   //    r = Sum_i rhoT q_i vec{g} . vec{W}_i,

   r = 0.0;
   if ( !traits_.dograv || !traits_.dofallout || traits_.dodry ) return;

   nhydro = qi.size(); // traits_.nlsector + traits_.nisector;

   tmp1 = utmp[GDIM];
   tmp2 = utmp[GDIM+1];

   // Set tmp space for gravity vector:
   for ( auto j=0; j<GDIM; j++ ) g[j] = utmp[j];
  
   // Compute Cartersian components of gravity vector field:
  *tmp1 = -GG;
   compute_vterm(*tmp1, g);

   r = 0.0;
   for ( auto j=0; j<nhydro; j++ ) {

     compute_vterm(*tvi[j], W_);       // compute W
     compute_fv(g, W_, *tmp1, *tmp2);  // g.W
     tmp2->pointProd(*qi[j]);          // qi * (g.W)
     r += *tmp2

   }
   r.pointProd(rhoT);                  // rhoT Sum_i q_i vec{g}.vec{W}_i

} // end of method compute_pe


//**********************************************************************************
//**********************************************************************************
// METHOD : assign_helpers
// DESC   : Assigns helper State components that are member data
//       
// ARGS   : u : incoming state vector
//          uf: incoming forcing state vector
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::assign_helpers(State &u, State &uf)
{
   GString     serr = "GMConv<TypePack>::assign_helpers: ";

   // Point helper arrays to corresp state components:
   qi_  = NULLPTR;
   tvi_ = NULLPTR;
   for ( auto j=0; j<nhydro_; j++ ) { //
     qi_ [j] = u[GDIM+3+j]; // set mass frac vector
     tvi_[j] = u[GDIM+3+nhydro_+j]; // set term speed for each qi
   }
   for ( auto j=0; j<GDIM; j++ ) fk_[j] = uf[j]; // kinetic forcing vector

   for ( auto j=0; j<GDIM            ; j++ ) s_    [j] = u[j];
   for ( auto j=0; j<traits_.nlsector; j++ ) qliq_ [j] = u[GDIM+3+j];
   for ( auto j=0; j<traits_.nisector; j++ ) qice_ [j] = u[GDIM+3+nliq+j];
   for ( auto j=0; j<traits_.nlsector; j++ ) tvliq_[j] = u[GDIM+3+  nliq+nice+j];
   for ( auto j=0; j<traits_.nisector; j++ ) tvice_[j] = u[GDIM+3+2*nliq+nice+j];

} // end of method assign_helpers


