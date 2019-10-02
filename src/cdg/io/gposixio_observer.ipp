//==================================================================================
// Module       : gposixio_observer.ipp
// Date         : 3/18/19 (DLR)
// Description  : Observer object for carrying out simple POSIX-based  
//                binary output.
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : ObserverBase.
//==================================================================================

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Instantiate with Traits
// ARGS   : traits: Traits sturcture
//**********************************************************************************
template<typename EquationType>
GPosixIOObserver<EquationType>::GPosixIOObserver(const EqnBasePtr &equation, Grid &grid,  typename ObserverBase<EquationType>::Traits &traits):
ObserverBase<EquationType>(equation, grid, traits),
bprgrid_         (TRUE),
bInit_          (FALSE),
cycle_              (0),
ocycle_             (0),
cycle_last_         (0),
time_last_        (0.0)
{ 
  this->traits_ = traits;
  this->grid_   = &grid;
} // end of constructor (1) method


//**********************************************************************************
//**********************************************************************************
// METHOD     : observe_impl
// DESCRIPTION: Prints state to files specified by traits. Format is:
//                  var1.CCCCCC.TTTTT.out,
//              where CCCCCC represents a cycle number, and TTTTT represents
//              the mpi task doing the writing.
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
void GPosixIOObserver<EquationType>::observe_impl(const Time &t, const State &u, const State &uf)
{
  init(t,u);

  mpixx::communicator comm;

  GIOTraits traits;
   
  if ( (this->traits_.itype == ObserverBase<EquationType>::OBS_CYCLE 
        && (cycle_-cycle_last_+1) >= this->traits_.cycle_interval)
    || (this->traits_.itype == ObserverBase<EquationType>::OBS_TIME  
        &&  t-time_last_ >= this->traits_.time_interval) 
    ||  cycle_ == 0 ) {
    traits.prgrid = bprgrid_;
    traits.wtime  = wtime_;
    traits.wtask  = wtask_;
    traits.wfile  = wfile_;
    traits.ivers  = ivers_;
    traits.dim    = GDIM;
    traits.gtype  = this->grid_->gtype();
    traits.index  = ocycle_;
    traits.cycle  = cycle_;
    traits.time   = t;
    traits.dir    = sodir_;
    gio_write_state(traits, *(this->grid_), u, state_index_, state_names_,  comm);
    gio_write_grid (traits, *(this->grid_), grid_names_,  comm);

    // Cycle through derived quantities, and write:
    print_derived(t, u, traits, comm);

    bprgrid_      = FALSE;
    cycle_last_   = cycle_;
    time_last_    = t;
    ocycle_++; // ouput cycle index
  }
  cycle_++;
  
} // end of method observe_impl


//**********************************************************************************
//**********************************************************************************
// METHOD     : init
// DESCRIPTION: Fill member index and name data based on traits
// ARGUMENTS  : t  : state time
//              u  : state variable
// RETURNS    : none.
//**********************************************************************************
template<typename EquationType>
void GPosixIOObserver<EquationType>::init(const Time &t, const State &u)
{

   if ( bInit_ ) return;

   char    stmp[1024];
   std::vector<GString> spref = {"x", "y", "z"};

   sidir_ = this->traits_.idir;
   sodir_ = this->traits_.odir;
   ivers_ = this->traits_.imisc;
   wtime_ = this->traits_.itag1;
   wtask_ = this->traits_.itag2;
   wfile_ = this->traits_.itag3;
 
   time_last_  = this->traits_.start_time ;
   ocycle_     = this->traits_.start_ocycle;
 
   // Set state index member data, if not already set:
   if ( state_index_.size()  <= 0 ) {
     if ( this->traits_.state_index.size() == 0 ) {
       for ( auto j=0; j<state_names_.size(); j++ ) {
         state_index_.push_back(j); 
       } 
     } 
     else {
       for ( auto j=0; j<this->traits_.state_index.size(); j++ ) {
         state_index_.push_back(this->traits_.state_index[j]); 
       } 
     }
   }

   // Set state names member data, if not already set:
   if ( state_names_.size()  <= 0 ) {
     if ( this->traits_.state_names.size() == 0 ) {
       for ( auto j=0; j<state_index_.size(); j++ ) {
         sprintf(stmp, "%s%d", "u", state_index_[j]+1);
         state_names_.push_back(stmp); 
       } 
     } 
     else {
       for ( auto j=0; j<state_index_.size(); j++ ) {
         state_names_.push_back(this->traits_.state_names[state_index_[j]].data()); 
       } 
     }
   }

   // Set grid names member data, if not already set:
   GINT ng = this->grid_->gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;
   if ( grid_names_.size()  <= 0 ) {
     if ( this->traits_.grid_names.size() == 0 ) {
       for ( auto j=0; j<ng; j++ ) {
         sprintf(stmp, "%sgrid", spref[j].c_str());
         grid_names_.push_back(stmp); 
       } 
     } 
     else {
       for ( auto j=0; j<ng; j++ ) {
         grid_names_.push_back(this->traits_.grid_names[j].data()); 
       } 
     }
   }

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
void GPosixIOObserver<EquationType>::print_derived(const Time &t, const State &u, GIOTraits &traits, const GC_COMM &comm)
{

  GINT               ntmp;
  GString            sop;   // math operation
  GTVector<GString>  sdqnames;
  GTVector<GINT>     iuin(3), iuout(3);
  State              tmp(3), uu(3), uout(3);
  char               stmp[1024];

    // Cycle through derived quantities, and write:
    for ( auto j=0; j<this->traits_.derived_quantities.size(); j++ ) {
      iuin    .resize(this->traits_.derived_quantities[j].icomponents.size());
      sdqnames.resize(this->traits_.derived_quantities[j].snames    .size());
      iuin     = this->traits_.derived_quantities[j].icomponents;
      sdqnames = this->traits_.derived_quantities[j].snames;
      sop      = this->traits_.derived_quantities[j].smath_op;
      uu.resize(iuin.size());
      ntmp     = this->utmp_->size() - uout.size();
      for ( auto i=0; i<uu  .size(); i++ ) uu  [i] = u[iuin[i]];
      for ( auto i=0; i<uout.size(); i++ ) uout[i] = (*(this->utmp_))[i];
      for ( auto i=0; i<uout.size(); i++ ) uout[i] = (*(this->utmp_))[i];
      for ( auto i=0; i<ntmp       ; i++ ) tmp [i] = (*(this->utmp_))[i+3];
   
      GMTK::domathop(*(this->grid_), uu, sop, tmp, uout, iuout);
      if ( sdqnames.size() < iuout.size() ) { // set default names
        sdqnames.resize(iuout.size());
        for ( auto i=0; i<iuout.size(); i++ ) {
          sprintf(stmp, "derived%d_c%d", j+1, i+1);
          sdqnames[j] = stmp; 
        }
      }

      gio_write_state(traits, *(this->grid_), uout, iuout, sdqnames,  comm);
    }


} // end of method print_derived

