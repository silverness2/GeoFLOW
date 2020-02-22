//==================================================================================
// Module       : gio_observer.ipp
// Date         : 3/18/19 (DLR)
// Description  : Observer object for carrying out binary output of state.
// Copyright    : Copyright 2019. Colorado State University. All rights reserved.
// Derived From : ObserverBase.
//==================================================================================

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Instantiate with EqnBasePtr, Grid, and Traits
// ARGS   : equation: EqnBasePtr
//          io_obj  : IO object
//          grid    : Grid object
//          traits  : Traits sturcture
//**********************************************************************************
template<typename EquationType>
GIOObserver<EquationType>::GIOObserver(const EqnBasePtr &equation, Grid &grid,  const IOBasePtr &io_ptr, typename ObserverBase<EquationType>::Traits &traits):
ObserverBase<EquationType>(equation, grid, traits),
bprgrid_         (TRUE),
bInit_          (FALSE),
cycle_              (0),
ocycle_             (0),
cycle_last_         (0),
time_last_        (0.0),
pIO_           (io_ptr)
{ 
  this->grid_  = &grid;
  stateinfo_   = equation->stateinfo(); 
//this->iotraits_  = &io_obj.get_traits(); 
} // end of constructor (1) method


//**********************************************************************************
//**********************************************************************************
// METHOD     : observe_impl
// DESCRIPTION: Prints state to files specified configured IO object.
//              NOTE: an internal cycle counter is maintained, as this 
//                    observer, like all others,  should be called at 
//                    each time step.
//
// ARGUMENTS  : t    : time, t^n, for state, uin=u^n
//              u    : state
//              uf   : forcing
//               
// RETURNS    : none.
//**********************************************************************************
template<typename EquationType>
void GIOObserver<EquationType>::observe_impl(const Time &t, const State &u, const State &uf)
{

  assert(bInit_ && "Object not initialized");

  mpixx::communicator comm;
  GINT                nstate=0;
  GTVector<GTVector<GFTYPE>>
                     *xnodes = &(this->grid_->xNodes());

  if ( (this->traits_.itype == ObserverBase<EquationType>::OBS_CYCLE 
        && (cycle_-cycle_last_+1) >= this->traits_.cycle_interval)
    || (this->traits_.itype == ObserverBase<EquationType>::OBS_TIME  
        &&  t-time_last_ >= this->traits_.time_interval) 
    ||  cycle_ == 0 ) {
    stateinfo_.sttype = 1; // 'state' state
    stateinfo_.nelems = this->grid_->nelems();
    stateinfo_.gtype  = this->grid_->gtype();
    stateinfo_.index  = ocycle_;
    stateinfo_.cycle  = cycle_;
    stateinfo_.time   = t;
    stateinfo_.svars  = this->traits_.state_names;

    for ( auto j=0; j<u.size(); j++ ) nstate += (stateinfo_.icomptype[j] != GSC_PRESCRIBED 
                                             &&  stateinfo_.icomptype[j] != GSC_NONE);
    up_.resize(nstate);
    for ( auto j=0; j<u.size(); j++ ) {
      if ( stateinfo_.icomptype[j] != GSC_PRESCRIBED
        && stateinfo_.icomptype[j] != GSC_NONE ) up_[j] = u[j];
    }
    pIO_->write_state(this->traits_.agg_state_name, stateinfo_, up_);

    if ( bprgrid_ ) {
      gridinfo_.sttype = 0; // grid 'state'
      gridinfo_.nelems = stateinfo_.nelems;
      gridinfo_.gtype  = stateinfo_.gtype;
      gridinfo_.index  = ocycle_;
      gridinfo_.cycle  = cycle_;
      gridinfo_.time   = t;
      gridinfo_.svars  = this->traits_.grid_names;
      gridinfo_.porder.resize(stateinfo_.porder.dim(1),stateinfo_.porder.dim(2)); 
      gridinfo_.porder = stateinfo_.porder;
      
      for ( auto j=0; j<gp_.size(); j++ ) gp_[j] = &(*xnodes)[j];
      pIO_->write_state(this->traits_.agg_grid_name, gridinfo_, gp_);
      bprgrid_ = FALSE;
    }

    // Cycle through derived quantities, and write:
    print_derived(t, u);

    cycle_last_   = cycle_;
    time_last_    = t;
    ocycle_++; // ouput cycle index
  }
  cycle_++;
  
} // end of method observe_impl


//**********************************************************************************
//**********************************************************************************
// METHOD     : init
// DESCRIPTION: Set member data based on state info
// ARGUMENTS  : info : state info
// RETURNS    : none.
//**********************************************************************************
template<typename EquationType>
void GIOObserver<EquationType>::init(StateInfo &info)
{

   char    stmp[1024];
   std::vector<GString> spref = {"x", "y", "z"};

   time_last_  = info.time ;
   ocycle_     = info.cycle;
 
   bInit_ = TRUE;

} // end of method init


//**********************************************************************************
//**********************************************************************************
// METHOD     : print_derived 
// DESCRIPTION: Write derived quantities
// ARGUMENTS  : t     : state time
//              u     : state variable used to compute derived quantities
//              traits: GIOTraits structure  for printing
//              comm  : communicator
// RETURNS    : none.
//**********************************************************************************
template<typename EquationType>
void GIOObserver<EquationType>::print_derived(const Time &t, const State &u)
{

  GINT               ntmp, nstate;
  GString            sop;   // math operation
  GTVector<GString>  sdqnames;
  GTVector<GINT>     iuin(3), iuout(3);
  State              tmp(3), uu(3), uout(3);
  GString            agg_derived;
  char               stmp[1024];

    // Cycle through derived quantities, and write:
    for ( auto j=0; j<this->traits_.derived_quantities.size(); j++ ) {
      iuin    .resize(this->traits_.derived_quantities[j].icomponents.size());
      sdqnames.resize(this->traits_.derived_quantities[j].snames    .size());
      iuin       = this->traits_.derived_quantities[j].icomponents;
      sdqnames   = this->traits_.derived_quantities[j].snames;
      agg_derived= this->traits_.derived_quantities[j].agg_sname;
      sop        = this->traits_.derived_quantities[j].smath_op;

      assert( iuin.size() > 0 && "Derived quantities require state component(s)");
      assert( iuin.min() >= 0 && iuin.max()< u.size()  && "Invalid component indices");
      if ( "" == sop ) continue; // nothing to do

      uu.resize(iuin.size());
      ntmp     = this->utmp_->size() - uout.size();
      for ( auto i=0; i<uu  .size(); i++ ) uu  [i] = u[iuin[i]];
      for ( auto i=0; i<uout.size(); i++ ) uout[i] = (*(this->utmp_))[i];
      for ( auto i=0; i<uout.size(); i++ ) uout[i] = (*(this->utmp_))[i];
      for ( auto i=0; i<3          ; i++ ) tmp [i] = (*(this->utmp_))[i+3];
   
      GMTK::domathop(*(this->grid_), uu, sop, tmp, uout, iuout);
      assert(sdqnames.size() >= iuout.size());
      for ( auto i=0; i<iuout.size(); i++ ) {
        this->grid_->get_ggfx().doOp(*uout[i], GGFX_OP_SMOOTH);
      }
      // Rest of stateinfo_ shoule have been set before call:
      stateinfo_.svars  = this->traits_.state_names;
      up_.resize(iuout.size());
      for ( auto j=0; j<up_.size(); j++ ) {
        up_[j] = uout[iuout[j]];
      }
      pIO_->write_state(agg_derived, stateinfo_, up_);
    }


} // end of method print_derived

