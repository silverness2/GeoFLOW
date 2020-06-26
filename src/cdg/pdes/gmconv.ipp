//==================================================================================
// Module       : gmconv.hpp
// Date         : 6/11/20 (DLR)
// Description  : Object defining a moist convection solver:
//
//                PDEs:
//                     TBD
//
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
#include <cstdlib>
#include <memory>
#include <cmath>
#include "gab.hpp"
#include "gext.hpp"
#include "gbdf.hpp"
#include "gexrk_stepper.hpp"

using namespace std;

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
bupdatebc_               (FALSE),
bsteptop_                (FALSE),
bvterm_                  (FALSE),
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
  GridIcos *icos = dynamic_cast<GGridIcos*>(grid_);

  assert(tmp.size() >= req_tmp_size() && "Insufficient tmp space provided");
  assert(!(GDIM==2 && icos!=NULLPTR && 2D spherical grid not allowed");

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
  if ( gmass_   != NULLPTR ) delete gmass_;
  if ( gimass_  != NULLPTR ) delete gimass_;
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

  GINT    ibeg, , nmfrac;
  GString serr = "GMConv<TypePack>::dudt_impl: ";

  // NOTE:
  // Make sure that, in init(), Helmholtz op is using only
  // weak Laplacian (q * mass term isn't being used), or there will 
  // be problems. This is required for explicit schemes, for
  // which this method is called.

  assert( !traits_.bsonserved ); // don't allow conservative form yet

  nmfrac = traits_.dodry ? 0 : traits_.nlsector + traits_.nisector + 1;

  // If non-conservative, compute RHS from:
  //     du/dt = -c(u).Grad u + nu nabla^2 u 
  // for each u, where c may be indep of u (pure advection):

  *urhstmp_[0] = *u[GDIM+1]; urhstmp[0]->rpow(-1.0) // 1/total mass

  // Mass fraction equations. These are:
  //   d rho_T + Div (rho_T v) = 0, 
  // for total mass
  //   d (dhot_T q_v) + Div (rho_T q_v v) = dot(s_v), 
  // for vapor mass fraction
  //   d (dhot_T q_i) + Div (rho_T q_i v) = dot(s_i) - , 
  // for vapor mass fraction
  // For non-conserved mass fractions, we solve
  //  df/dt u.Grad f = dot(s)/rho_tot
  ibeg   = GDIM + 2;
  for ( auto k=ibeg; k<ibeg+nmfrac; k++ ) {
    // Multiply mass frac by total den:
    *urhstmp_[1] = (*u[k]) * (*urhstmp[0]); // rho_T * q_i
    
  }
  
  // Energy equation:
  // Momentum equations:
  for ( auto k=0; k<GDIM; k++ ) {
    gadvect_->apply(*u[k], u, uoptmp_, *dudt[k]);     // apply advection
    ghelm_->opVec_prod(*u[k], uoptmp_, *urhstmp_[0]); // apply diffusion
    GMTK::saxpby<GFTYPE>(*urhstmp_[0], -1.0, *dudt[k], -1.0);
    gimass_->opVec_prod(*urhstmp_[0], uoptmp_, *dudt[k]); // apply M^-1
    if ( traits_.bforced && uf[k] != NULLPTR ) *dudt[k] += *uf[k];
  }
  
} // end of method dudt_impl


//**********************************************************************************
//**********************************************************************************
// METHOD : step_impl (1)
// DESC   : Step implementation method entry point
// ARGS   : t   : time
//          uin : input state, modified on output with update
//                If doing pure advection, this state contains the *unevolved*
//                advection velocities, ci, which must be separated from the 
//                single evolved state variable:
//                     uin = [u_evolve, c1, c2, c3]
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
//                If doing pure advection, this state contains the *unevolved*
//                advection velocities, which must be separated from the 
//                single evolved state variable.
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

  // Set dissipation from traits. Note that
  // this is set to be constant, based on configuration,
  // even though the solver can accommodate spatially
  // variable dissipation:
  nu_.resize(1);
  nu_ = traits_.nu;

  // Find no. state and solve members, and component types:
  for( auto j=0; j<GDIM; j++ ) icomptype->push_back(GSC_KINETIC); 
  icomptype->push_back(GSC_TOTDENSITY); 
  icomptype->push_back(GSC_INTENERGY); 
  if ( traits_.dodry  {
    traits_.nsolve = GDIM + 2;
    traits_.nstate = traits_.nsolve;
  }
  else {
    n = traits_.nlsector + traits_.nisector;
    traits_.nsolve = GDIM + 2 + n;
    traits_.nstate = traits_.nsolve;
    for( auto j=0; j<traits_.nlsector; j++ ) icomptype->push_back(GSC_LIQDENSITY); 
    for( auto j=0; j<traits_.nisector; j++ ) icomptype->push_back(GSC_ICEDENSITY); 
    if traits_.dofallout ) { // include terminal velocities:
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
  gmass_   = new GMass(*grid_);
  ghelm_   = new GHelmholtz(*grid_);

  ghelm_->set_Lap_scalar(nu_);

  
  if ( traits_.isteptype ==  GSTEPPER_EXRK ) {
    gimass_ = new GMass(*grid_, TRUE); // create inverse of mass
  }

  // If doing semi-implicit time stepping; handle viscous term 
  // (linear) inplicitly, which implies using full Helmholtz operator:
  if ( traits_.isteptype == GSTEPPER_BDFAB || traits_.isteptype == GSTEPPER_BDFEXT ) {
    assert(FALSE && "Implicit time stepping not yet supported");
  }

  if ( traits_.bconserved && !doheat_ ) {
    assert(FALSE && "Conservation not yet supported");
    gpdv_  = new GpdV(*grid_,*gmass_);
//  gflux_ = new GFlux(*grid_);
    assert( (gmass_   != NULLPTR
          && ghelm_   != NULLPTR
          && gpdv_    != NULLPTR) && "1 or more operators undefined");
  }
  if ( !traits_.bconserved && !doheat_ ) {
    gadvect_ = new GAdvect(*grid_);
    assert( (gmass_   != NULLPTR
          && ghelm_   != NULLPTR
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
  ggrid_->minlength(&dxmin_);

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
void GMConv<TypePack>::apply_bc_impl(const Time &t, State &u, const State &ub)
{
  GTVector<GTVector<GSIZET>>  *igbdy = &grid_->igbdy_binned();
  GTVector<GTVector<GSIZET>>  *ilbdy = &grid_->ilbdy_binned();

  // Use indirection to set the global field node values
  // with domain boundary data. ub must be updated outside 
  // of this method.

  // NOTE: This is useful to set Dirichlet-type bcs only. 
  // Neumann bcs type have to be set with the
  // differential operators themselves, though the node
  // points in the operators may still be set from ub
 
  GSIZET   ib;
  GBdyType itype; 
  for ( GSIZET m=0; m<igbdy->size(); m++ ) { // for each type of bdy in gtypes.h
    itype = static_cast<GBdyType>(m);
    if (// itype == GBDY_NEUMANN
         itype == GBDY_PERIODIC
     ||  itype == GBDY_OUTFLOW
     ||  itype == GBDY_SPONGE 
     ||  itype == GBDY_NONE   ) continue;
    for ( GSIZET k=0; k<u.size(); k++ ) { // for each state component
      if ( ub[k] == NULLPTR ) continue;
      for ( GSIZET j=0; j<(*igbdy)[m].size(); j++ ) { // set Dirichlet-like value
        ib = (*igbdy)[m][j];
        (*u[k])[ib] = (*ub[k])[j];
      } 
    } 
  } 

  // Handle 0-Flux bdy conditions. This
  // is computed by solving
  //    vec{n} \cdot vec{u} = 0
  // for 'dependent' component set in grid.
  //
  // Note: We may want to switch the order of the
  //       following loops to have a better chance
  //       of vectorization. Unrolling likely
  //       won't occur:
  GINT                         id;
  GSIZET                       il;
  GFTYPE                       sum, xn;
  GTVector<GTVector<GFTYPE>>  *n    = &grid_->bdyNormals();
  GTVector<GINT>              *idep = &grid_->idepComp  ();

  itype = GBDY_0FLUX;
  for ( auto j=0; j<(*igbdy)[itype].size(); j++ ) { 
    ib = (*igbdy)[itype][j]; // index into vector array
    il = (*ilbdy)[itype][j]; // index into bdy array (for normals, e.g.)
    id = (*idep)[ib];        // dependent vector component
    xn = (*n)[id][ib];       // n_id == normal component for dependent vector comp
    sum = 0.0;
    for ( auto k=0; k<u.size(); k++ ) { // for each vector component
      if ( k != (*idep)[ib] ) sum -= (*n)[k][il] * (*u[k])[ib];
    }
    (*u[id])[ib] = sum / xn;
  }

  // Handle no-slip bdy conditions.
  itype = GBDY_NOSLIP;
  for ( auto k=0; k<u.size(); k++ ) { // for each velocity component
    for ( auto j=0; j<(*igbdy)[itype].size(); j++ ) { 
      ib = (*igbdy)[itype][j]; // index into vector array
      (*u[k])[ib] = 0.0;
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
//          utmp : tmp vectors; at least 1 required; only first one used
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
//          utmp : tmp vectors; at least 2 required; only first two used
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
//          utmp : tmp vectors; at least 3 required; only first 3 used
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
// METHOD : compute_div
// DESC   : Compute flux divergence  of quantity q:
//             Div ( q v)
//          with the intent of handling this conservatively or not.  
// ARGS   : q    : quantiy whose flux we compute divergence of
//          v    : vector of velocity components
//          utmp : tmp vectors; at least GDIM+1 required; only first GDIM+1 used
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
// DESC   : Compute 'terminal' velocity state vector velocity components  from 
//          prescribed terminal velocity data and store for use in W_ member array
// ARGS   : u     : total state vector
//          ihydro: index of hydrometeor (0 to nliquids + nices - 1)
//          utmp  : tmp vector array
// RETURNS: none. Member data, W_, is set here
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::compute_vterm(State &u, GINT ihydro, State &utmp)
{
   GString    serr = "GMConv<TypePack>::compute_vterm: ";
   GINT       ibeg = GDIM+3, nhydro;
   Value      r, x, y, z;
   Value      lat, long;
   GTVector<GTVector<GFTYPE>> 
             *xnodes = &grid.xNodes();
   StateComp *tvi; // ihydro-th terminal velocity

   GGridIcos *icos = dynamic_cast<GGridIcos*>(grid_);
   GGridIcos *box  = dynamic_cast <GGridBox*>(grid_);

   nhydro = traits_.nlsector + traits_.nisector;
   if ( nhydro == 0 ) return; // return if there are no hydrometeors

   assert( ihydro >= 0 && ihydro < nhydro );

   // Check if already computed:
   if ( traits_.bvarvterm && bvterm_ ) return;


   // Specify components for Cartesian grid:
   if ( box != NULLPTR ) {
     // In Cartesian coords, select the 'z' direction
     // as preferred 'fallout' direction. In 2d, this
     // will be the 2-coord; in 3d, the 3-coord:
     W_[GDIM-1] = &tvi; 
     bvterm_ = TRUE; // term vel. computed
     return;
   }

   // Specify components for spherical grid. Here,
   // the preferred direction is along radial, and
   // we assume that the prescribed term velocity
   // is the radial component only, which must be
   // transformed into Cartesian components:

   
   tvi = u[ibeg+ihydro];
   for ( auto i=0; i<(*xnodes)[0].size(); i++ ) {
      x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
      r     = sqrt(x*x + y*y + z*z);
      lat   = asin(z/r);
      long  = atan2(y,x);

     *W_[0]  = cos(lat)*cos(long);
     *W_[1]  = cos(lat)*sin(long);
     *W_[2]  = sin(lat);
   }

   for ( auto j=0; j<W_.size(); j++ ) {
     if ( tvi->size() == 1 ) {
       *W_[j] *= (*tvi)[0];
     }
     else {
       *W_[j] *= (*tvi);
     } 
   }

   bvterm_ = TRUE;

} // end of method compute_vterm


//**********************************************************************************
//**********************************************************************************
// METHOD : compute_falloutsrc
// DESC   : Compute total effect due to fallout of density quantity, g. 
//            r = Sum_i Div [q_i g vector(W_i)] 
//          where vector(W_i) is the fallout terminal velocity for precipitating
//          species i (from either liquid or ice sectors of state). 
// ARGS   : 
//          g    : quantiy to fall out ('flux-out'). May be: total density,
//                 momentum density component, internal energy density
//          qi   : hydrometeor mass fractions
//          tvi  : terminal velocity vector for each qi hydrometeor
//          jexcl: which hydrometeor index to exclude from sum. If jexcl<0,
//                 exclude none. This index must be the 0-starting index of
//                 the hydrometeor in the qi array.
//          utmp : tmp vectors; at least 2*GDIM+3 required
//          r    : fallout src field
// RETURNS: none.
//**********************************************************************************
template<typename TypePack>
void GMConv<TypePack>::compute_falloutsrc(StateComp &g, State &qi, State &tvi, GINT jexcl, State &utmp, StateComp &r)
{
   GString     serr = "GMConv<TypePack>::compute_falloutsrc: ";
   GINT        nhydro; // no. hydrometeors
   StateComp  *div, *qg; 
   State       vterm(GDIM); // term velocity vector

   assert(tvi.size() == qi.size());
   assert(utmp.size() >= 2*GDIM+1);
   assert(jexcl < 0 || (jexcl >=0 && jexcl < qi.size()));

   // Compute:
   //    r = -Sum_i Div (rho_t q_i vector(W)_i )

   r = 0.0;
   if ( !traits_.dofallout || traits_.dodry ) return;

   nhydro = traits_.nlsector + traits_.nisector;

   vterm = NULLPTR;
   for ( auto j=0; j<GDIM; j++ ) vterm[j] = utmp[GDIM+3+j];

   qg    = utmp[GDIM+1];    // temp
   div   = utmp[GDIM+2];    // temp
   for ( auto j=0; j<nhydro; j++ ) {
     if ( j == jexcl ) continue;  

     *qg = g; (*qg) *= (*qi[j]); // compute g q_i

     // Convert terminal velocities to required 
     // (Cartesian) components. vterm may point
     // to tvi, or to v_ member data, or may be NULL:
     compute_vterm(*tvi[j], utmp, vterm);

     // Compute i_th contribution to source term:
     compute_div(*qg, vterm, utmp, *div);
     r += *div;
   }

} // end of method compute_falloutsrc


