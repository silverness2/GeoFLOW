//==================================================================================
// Module       : gupdatebdy_user.cpp
// Date         : 7/11/19 (DLR)
// Description  : Boundary update function implementations specified by
//                user
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#include "gupdatebdy_user.hpp"


namespace gupdatebdy {



//**********************************************************************************
//**********************************************************************************
// METHOD : impl_simple_outflow
// DESC   : Update bdy using simple outflow. Applies
//          only to GBDY_OUTFLOW bdy types. Called each time step.
// ARGS   : ptree : main prop tree
//          sconfiig : config block name
//          grid  : grid
//          stinfo: StateInfo
//          t     : time
//          utmp  : tmp arrays; require at least ub.size arrays
//          u     : current state, unchanged here
//          ub    : bdy vectors (one for each state element, unless NULL)
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_simple_outflow(const PropertyTree &ptree, GString &sconfig, GGrid &grid, StateInfo &stinfo, Time &time, State &utmp, State &u, State &ub)
{
  Time             tt = time;
  State            uu(u.size());
  GString          serr = "impl_simple_outflow: ";

  GTVector<GTVector<GSIZET>> *igbdy = &grid.igbdy_binned();

  // Set from State vector, u:
  for ( auto k=0; k<u.size(); k++ ) {
    for ( auto j=0; j<(*igbdy)[GBDY_OUTFLOW].size()
       && ub[k] != NULLPTR; j++ ) {
      (*ub[k])[j] = (*u[k])[(*igbdy)[GBDY_OUTFLOW][j]];
    }
  }

  return TRUE;

} // end, impl_simple_outflow


//**********************************************************************************
//**********************************************************************************
// METHOD : impl_sphere_sponge
// DESC   : Update a outer spherical sponge layer.  
//          Unlike for other bdy conditions, this one, if called
//          will be applied for the entire computational volume, 
//          focusing on flow in the far field. In other words, 
//          it is appled to all points with r > rs below, and
//          does not care about what the spherical outer boundary 
//          conditions are set to.  Called each time step.
// ARGS   : ptree    : main prop tree
//          sconfiig : config block name
//          grid     : grid
//          stinfo   : StateInfo
//          t        : time
//          utmp     : tmp arrays; require at least ub.size arrays
//          u        : current state, unchanged here
//          ub       : bdy vectors (one for each state element, unless NULL)
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_sphere_sponge(const PropertyTree &ptree, GString &sconfig, GGrid &grid, StateInfo &stinfo, Time &time, State &utmp, State &u, State &ub)
{
  Time             tt = time;
  GINT             idstate;
  GFTYPE           beta, exponent, ifact, ro, rs, sig0;
  GFTYPE           r, x, y, z;
  State            uu(u.size());
  std::vector<GINT>
                   istate;
  std::vector<GFTYPE> 
                   farfield;
//GTVector<GTVector<GINT>> 
//                *igbdycf = &grid.igbdycf_binned(); 
//GTVector<GTVector<GSIZET>> 
//                *igbdy = &grid.igbdy_binned();

  GTVector<GTVector<GFTYPE>> 
                  *xnodes = &grid.xNodes();
  PropertyTree     sptree; // sponge config tree
  GString          serr = "impl_sphere_sponge: ";

  assert(GDIM == 3);

  // Get parameters from ptree:
  sptree   = ptree.getPropertyTree(sconfig);
  ro       = sptree.getValue<GFTYPE>("radiuso");       // outer grid radius
  rs       = sptree.getValue<GFTYPE>("layer_radius");  // inner sponge radius
  exponent = sptree.getValue<GFTYPE>("exponent", 4);   // exponent 
  sig0     = sptree.getValue<GFTYPE>("sigma", 1.0);    // diffusion constant
  istate   = sptree.getArray  <GINT>("state_index");   // state ids to 'sponge'
  farfield = sptree.getArray<GFTYPE>("far_field");     // far-field state values

  ifact    = 1.0/(ro - rs);

  // This method applies a sponge layer to only the outer
  // part of a spherical grid

  // Update state due to sponge layer:
  // Note: This is equiavalent to adding a dissipation 
  //       term to the RH of the operator-split equation, s.t.:
  //        du/dt = -sig(r) (u - u_infinity)
  //       where
  //        sig(r) = sig_0 [(r - rs)/(ro - rs)]^exponent
  //       and u_infinity is the far-field solution
  // Note: We may have to re-form this scheme if we use semi-implicit
  //       or implicit time stepping methods!
  for ( auto k=0; k<istate.size(); k++ ) { // for each state component
    idstate = istate[k];
    for ( auto j=0; j<u[idstate]->size(); j++ ) {
      x    = (*xnodes)[0][k]; y = (*xnodes)[1][k]; z = (*xnodes)[2][k];
      r    = sqrt(x*x + y*y + z*z); 
//    igb  = (*igbdy)  [GBDY_SPONGE][j];
//    igf  = (*igbdycf)[GBDY_SPONGE][j];
      beta = r >= rs ? pow(ifact*(r-rs),exponent) : 0.0; // check if in sponge layer
//    (*u[idstate])[j]  = (1.0-beta)*(*u[idstate])[j] + beta*farfield[k];
      (*u[idstate])[j] -= sig0*beta*( (*u[idstate])[j] - farfield[k] );
    }
  }

  return TRUE;

} // end, impl_sphere_sponge


} // end, gupdatebdy namespace
