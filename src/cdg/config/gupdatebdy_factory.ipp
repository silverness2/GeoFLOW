//==================================================================================
// Module       : gupdatebdy_factory
// Date         : 7/11/19 (DLR)
// Description  : GeoFLOW state initialization factory
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================

//**********************************************************************************
//**********************************************************************************
// METHOD : update
// DESC   : Do bdy update
// ARGS   : ptree  : main property tree
//          grid   : GGrid operator
//          stinfo : StateInfo
//          time   : initialization time
//          utmp   : tmp arrays
//          u      : state to be initialized. 
//          ub     : boundary state 
// RETURNS: none.
//**********************************************************************************
template<typename EquationType>
GBOOL GUpdateBdyFactory<EquationType>::update(const PropertyTree& ptree, Grid &grid, StateInfo &stinfo, Time &time, State &utmp, State &u, State &ub)
{
  GBOOL         bret = FALSE, use_inits;
  State         uu(u.size());
  GString       sgrid, supdate;
  PropertyTree  gtree;
  std::function<void(const geoflow::tbox::PropertyTree& ptree,GString &supdate, Grid &grid,
                     StateInfo &stinfo, Time &time, State &utmp, State &u, State &ub)>
                mycallback;
  sgrid    = ptree.getValue<GString>("grid_type");
  gtree    = ptree.getPropertyTree(sgrid);
  supdate  = gtree.getValue<GString>("bdy_update_method","none");
  use_inits= gtree.getValue<GBOOL>  ("use_state_init",FALSE);

  if ( "updateb_none" == supdate
    || "none"         == supdate
    || ""             == supdate ) {
    bret = TRUE;
  }
  else if ( use_inits ) {
    assert(utmp.size() >= 2*u.size() && "Tmp array is too small!");
    for ( auto i=0; i<uu.size(); i++ ) {
      uu[i] = utmp[u.size()+i];
     *uu[i] = *u[i];
    }
    if ( bret ) {
       mycallback =
                    [](
                     const geoflow::tbox::PropertyTree& ptree,GString &supdate, Grid &grid, 
                     StateInfo &stinfo, Time &time, State &utmp, State &u, State &ub)
                     {set_bdy_from_state(ptree, supdate, grid, stinfo, time, utmp, u, ub);};
      grid.set_update_bdy_callback(mycallback);
    }

  }
  else if ( "simple_outflow" == supdate ) {
//  bret = gupdatebdy::impl_simple_outflow (ptree, grid, time, utmp, u, ub);
    mycallback =
                  [](
                  const geoflow::tbox::PropertyTree& ptree, GString &supdate, Grid &grid, 
                  StateInfo &stinfo, Time &time, State &utmp, State &u, State &ub)
                  {gupdatebdy::impl_simple_outflow(ptree, supdate, grid, stinfo, time, utmp, u, ub);};
    grid.set_update_bdy_callback(mycallback);
  }
  else if ( "sphere_sponge" == supdate ) {
    mycallback
         = [](
           const geoflow::tbox::PropertyTree& ptree, GString &supdate, Grid &grid, 
           StateInfo &stinfo, Time &time, State &utmp, State &u, State &ub)
           {gupdatebdy::impl_sphere_sponge(ptree, supdate, grid, stinfo, time, utmp, u, ub);};
    grid.set_update_bdy_callback(mycallback);
  }
  else {
    assert(FALSE && "Specified bdy update method unknown");
  }

  return bret;

} // end, init method update


//**********************************************************************************
//**********************************************************************************
// METHOD : set_bdy_from_state
// DESC   : use state var, u, to set bdy, ub
// ARGS   : ptree  : main property tree
//          sconfig: config block name
//          grid   : GGrid operator
//          stinfo : StateInfo
//          time   : initialization time
//          utmp   : tmp arrays
//          u      : state to be initialized. 
//          ub     : boundary state 
// RETURNS: none.
//**********************************************************************************
template<typename EquationType>
void GUpdateBdyFactory<EquationType>::set_bdy_from_state(const geoflow::tbox::PropertyTree& ptree, GString &sconfig, Grid &grid, StateInfo &stinfo, Time &time, State &utmp, State &u, State &ub)
{
  GBOOL         bret=FALSE;
  GBOOL         use_inits; // use state init method to set bdy?
  GString       sgrid, supdate;

  GTVector<GTVector<GSIZET>> *igbdy = &grid.igbdy_binned();

  bret = GInitStateFactory<EquationType>::init(ptree, grid, stinfo, time, utmp, ub, u);

  // Set from State vector, u and others that we _can_ set:
  for ( auto k=0; k<u.size(); k++ ) {
    for ( auto j=0; j<(*igbdy)[GBDY_DIRICHLET].size()
       && ub[k] != NULLPTR; j++ ) {
      (*ub[k])[j] = (*u[k])[(*igbdy)[GBDY_DIRICHLET][j]];
    }
    for ( auto j=0; j<(*igbdy)[GBDY_INFLOWT].size()
       && ub[k] != NULLPTR; j++ ) {
      (*ub[k])[j] = (*u[k])[(*igbdy)[GBDY_INFLOWT][j]];
    }
  }

} // end, set_bdy_from_state

