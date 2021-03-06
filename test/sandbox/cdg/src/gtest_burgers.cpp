//==================================================================================
// Module       : gtest_burgers.cpp
// Date         : 12/24/18 (DLR)
// Description  : GeoFLOW test of Burgers PDE
// Copyright    : Copyright 2018-19. Colorado State University. All rights reserved.
// Derived From : 
//==================================================================================
#include "gexec.h"
#include "gtypes.h"
#include <cstdio>
#include <cmath>
#include <unistd.h>
#include <iostream>
#include <memory>
#include <cstdlib>
#include <cassert>
#include <random>
#include <typeinfo>
#include "gcomm.hpp"
#include "ggfx.hpp"
#include "gllbasis.hpp"
#include "gmorton_keygen.hpp"
#include "gburgers.hpp"
#include "ggrid_box.hpp"
#include "ggrid_factory.hpp"
#include "gposixio_observer.hpp"
#include "pdeint/equation_base.hpp"
#include "pdeint/integrator_factory.hpp"
#include "pdeint/mixer_factory.hpp"
#include "pdeint/observer_factory.hpp"
#include "pdeint/null_observer.hpp"
#include "tbox/property_tree.hpp"
#include "tbox/mpixx.hpp"
#include "tbox/global_manager.hpp"
#include "tbox/input_manager.hpp"
#include "gtools.h"
//#include "gio.h"

#if defined(_G_USE_GPTL)
    #include "gptl.h"
#endif

using namespace geoflow::pdeint;
using namespace geoflow::tbox;
using namespace std;


template< // default template arg types
typename StateType = GTVector<GTVector<GFTYPE>*>,
typename GridType  = GGrid,
typename ValueType = GFTYPE,
typename DerivType = StateType,
typename TimeType  = ValueType,
typename JacoType  = StateType,
typename SizeType  = GSIZET,
>
struct EquationTypes {
        using State      = StateType;
        using Grid       = GridType;
        using Value      = ValueType;
        using Derivative = DerivType;
        using Time       = TimeType;
        using Jacobian   = JacoType;
        using Size       = SizeType;
};

GBOOL  bench_=FALSE;
GINT   nstate_=GDIM;  // number 'state' arrays
GINT   nsolve_=GDIM;  // number *solved* 'state' arrays
GGrid *grid_;
GTVector<GTVector<GFTYPE>*> u_;
GTVector<GTVector<GFTYPE>*> c_;
GTVector<GTVector<GFTYPE>*> ua_;
GTVector<GTVector<GFTYPE>*> ub_;
GTVector<GTVector<GFTYPE>*> uf_;
GTVector<GTVector<GFTYPE>*> utmp_;
GTVector<GTVector<GFTYPE>*> uold_;
GTVector<GFTYPE> nu_(3);
PropertyTree ptree;       // main prop tree
GC_COMM      comm_ ;      // communicator


using MyTypes = EquationTypes<>;       // Define types used
using EqnBase = EquationBase<MyTypes>; // Equation base type
using EqnImpl = GBurgers<MyTypes>;     // Equation implementation
using MixBase= MixerBase<MyTypes>;  // Mixing base Type
using MixBasePtr 
              = std::shared_ptr<MixBase>; // Mixing base ptr
using ObsBase = ObserverBase<EqnBase>; // Observer base Type

void compute_analytic(GGrid &grid, GFTYPE &t, const PropertyTree& ptree, GTVector<GTVector<GFTYPE>*> &u);
void update_dirichlet(const GFTYPE &t, GTVector<GTVector<GFTYPE>*> &u, GTVector<GTVector<GFTYPE>*> &ub);
void steptop_callback(const GFTYPE &t, GTVector<GTVector<GFTYPE>*> &u, const GFTYPE &dt);
void create_observers(PropertyTree &ptree, GSIZET icycle, GFTYPE time, 
std::shared_ptr<std::vector<std::shared_ptr<ObserverBase<MyTypes>>>> &pObservers);
void create_mixer(PropertyTree &ptree, MixBasePtr &pMixer);
void gresetart(PropertyTree &ptree);
void do_bench(GString sbench, GSIZET ncyc);

//#include "init_pde.h"


int main(int argc, char **argv)
{

    GString serr ="main: ";
    GBOOL  bret, dopr=TRUE;
    GBOOL  bgridwritten=FALSE;
    GINT   iopt;
    GINT   ilevel=0;// 2d ICOS refinement level
    GINT   np=1;    // elem 'order'
    GSIZET itindex=0;    // restart flag/index
    GSIZET icycle=0;       // curr time cycle
    GFTYPE tt;
    std::vector<GINT> ne(3); // # elements in each direction in 3d
    std::vector<GINT> pstd(GDIM);  // order in each direction
    GString sgrid;// name of JSON grid object to use
    GTVector<GString> savars, scvars, sdvars;
    GTMatrix<GINT> p; // needed for restart, but is dummy
    char stmp[1024];

    typename MyTypes::Time  t  = 0;
    typename MyTypes::Time  dt = 0.1;

    // Initialize comm & global environment:
    mpixx::environment env(argc,argv); // init GeoFLOW comm
    mpixx::communicator world;
    GINT                myrank;
    GINT                ntasks;
    GlobalManager::initialize(argc,argv); 
    GlobalManager::startup();
    comm_ = world; // need this for solver(s) & grid
    myrank = world.rank();
    ntasks = world.size();


    // Read main prop tree; may ovewrite with
    // certain command line args:
    PropertyTree eqptree;     // equation props
    PropertyTree gridptree;   // grid props
    PropertyTree stepptree;   // stepper props
    PropertyTree dissptree;   // dissipation props
    PropertyTree tintptree;   // time integration props

    EH_MESSAGE("main::call load_file...");
    ptree    = InputManager::getInputPropertyTree();       // main param file structure
    itindex = ptree.getValue<GSIZET>   ("restart_index");

    EH_MESSAGE("main: load_file done.");
    // Create other prop trees for various objects:
    sgrid       = ptree.getValue<GString>("grid_type");
    pstd        = ptree.getArray<GINT>("exp_order");
    eqptree     = ptree.getPropertyTree("burgers_traits");
    gridptree   = ptree.getPropertyTree(sgrid);
    stepptree   = ptree.getPropertyTree("stepper_props");
    dissptree   = ptree.getPropertyTree("dissipation_traits");
    tintptree   = ptree.getPropertyTree("time_integration");

    ne          = gridptree.getArray<GINT>("num_elems");  // may be modified by command line
#if 1

    // Parse command line. ':' after char
    // option indicates that it takes an argument.
    // Note: -i reserved for InputManager:
    while ((iopt = getopt(argc, argv, "i:j:k:l:m:np:h")) != -1) {
      switch (iopt) {
      case 'i': // handled by InputManager
          break;
      case 'j': // get # elements in r/x
          ne[0] = atoi(optarg);
          gridptree.setArray<GINT>("num_elems",ne);
          break;
      case 'k': // get # elements in lat/y
          ne[1] = atoi(optarg);
          gridptree.setArray<GINT>("num_elems",ne);
          break;
      case 'm': // get # elements in long/z
          ne[2] = atoi(optarg);
          gridptree.setArray<GINT>("num_elems",ne);
          break;
      case 'l': // # 2d refinement level
          ilevel = atoi(optarg);
          gridptree.setValue<GINT>("ilevel",ilevel);
          break;
      case 'n': // no output
          dopr = FALSE;
          break;
      case 'p': // get nodal exp order
          np = atoi(optarg);
          pstd.assign(GDIM, np);
          break;
      case 'h': // help
          std::cout << "usage: " << std::endl <<
          argv[0] << " [-h] [-j #Elems in x/r] [-k #Elems in y/lat]  [-m #Elems in z/long] [-l refine level] -p expansion order] " << std::endl;
          std::cout << "Note: Icos grid only requires refine level to specify number of elements. " << std::endl;
          exit(1); 
          break;
      case ':': // missing option argument
          std::cout << argv[0] << ": option " << optopt << " requires an argument" << std::endl;
          exit(1); 
          break;
      case '?':
      default: // invalid option
          std::cout << argv[0] << ": option " << optopt << " invalid" << std::endl;
          exit(1);
          break;
      }
    }

#endif
    ptree.setValue<GString>("exp_order_type","constant");
    bench_ = ptree.getValue<GBOOL>("benchmark");

    // Set solver traits from prop tree:
    GFTYPE nu_scalar;
    GBurgers<MyTypes>::Traits solver_traits;
    solver_traits.doheat     = eqptree  .getValue<GBOOL>  ("doheat", FALSE);
    solver_traits.bpureadv   = eqptree  .getValue<GBOOL>  ("bpureadv", FALSE);
    solver_traits.bconserved = eqptree  .getValue<GBOOL>  ("bconserved", FALSE);
    solver_traits.bforced    = eqptree  .getValue<GBOOL>  ("use_forcing", FALSE);
    solver_traits.variabledt = eqptree  .getValue<GBOOL>  ("variable_dt", FALSE);
    solver_traits.courant    = eqptree  .getValue<GFTYPE> ("courant", 0.5);
    solver_traits.ssteptype  = stepptree.getValue<GString>("stepping_method");
    solver_traits.itorder    = stepptree.getValue <GINT>  ("time_deriv_order");
    solver_traits.inorder    = stepptree.getValue <GINT>  ("extrap_order");
    nu_scalar                = dissptree.getValue<GFTYPE> ("nu");
    nu_.resize(1); 
    nu_ = nu_scalar; 
    
#if defined(_G_USE_GPTL)
    // Set GTPL options:
    GPTLsetoption (GPTLcpu, 1);
#endif
    // Initialize timer:
    GTimerInit();

    // Create basis:
    GTVector<GNBasis<GCTYPE,GFTYPE>*> gbasis(GDIM);
    for ( GSIZET j=0; j<GDIM; j++ ) {
      gbasis [j] = new GLLBasis<GCTYPE,GFTYPE>(pstd[j]);
    }
    
    EH_MESSAGE("main: build grid...");

    GTimerStart("gen_grid");

    // Create grid:
    grid_ = GGridFactory::build(ptree, gbasis, comm_);
    GComm::Synch(comm_);

    GTimerStop("gen_grid");


    EH_MESSAGE("main: initialize gather/scatter...");

    GTimerStart("init_ggfx_op");

    // Initialize gather/scatter operator:
    GGFX<GFTYPE> ggfx;
    init_ggfx(ptree, *grid_, ggfx);

    GTimerStop("init_ggfx_op");
    EH_MESSAGE("main: gather/scatter initialized.");

    // Create state and tmp space:
    EH_MESSAGE("main: set up tmp space...");
    if      ( solver_traits.doheat   
           || solver_traits.bpureadv ) { 
      nstate_ = sgrid == "grid_icos" ? 4 : GDIM + 1; nsolve_ = 1; // 1-state + no. v-components
    }
    utmp_.resize(24);
    uold_.resize(nsolve_);
    u_ .resize(nstate_);
    ua_.resize(nsolve_);
    ub_.resize(nsolve_);
    uf_.resize(nsolve_); uf_ = NULLPTR;
    c_ .resize(sgrid == "grid_icos" ? 3 : GDIM);
    for ( GSIZET j=0; j<uold_.size(); j++ ) uold_[j] = new GTVector<GFTYPE>(grid_->size());
    for ( GSIZET j=0; j<utmp_.size(); j++ ) utmp_[j] = new GTVector<GFTYPE>(grid_->size());
    for ( GSIZET j=0; j<u_   .size(); j++ ) u_   [j] = new GTVector<GFTYPE>(grid_->size());
    for ( GSIZET j=0; j<ua_  .size(); j++ ) ua_  [j] = new GTVector<GFTYPE>(grid_->size());
    for ( GSIZET j=0; j<ub_  .size(); j++ ) ub_  [j] = new GTVector<GFTYPE>(grid_->nbdydof());
    if ( solver_traits.bforced )
    for ( GSIZET j=0; j<ub_  .size(); j++ ) uf_  [j] = new GTVector<GFTYPE>(grid_->nbdydof());
    if ( solver_traits.bpureadv || solver_traits.doheat ) 
    for ( GSIZET j=0; j<c_   .size(); j++ ) c_   [j] = u_[j+1];

    // Create equation set:
    EH_MESSAGE("main: create equation...");
    std::shared_ptr<EqnImpl> eqn_impl(new EqnImpl(ggfx, *grid_, u_, solver_traits, utmp_));
    std::shared_ptr<EqnBase> eqn_base = eqn_impl;

    // Set Dirichlet bdy state update function:
    EH_MESSAGE("main: set bdy update functions...");
    std::function<void(const GFTYPE &t, GTVector<GTVector<GFTYPE>*> &u, 
                                        GTVector<GTVector<GFTYPE>*> &ub)>  
        fcallback = update_dirichlet; // set tmp function with proper signature for...
    std::function<void(const GFTYPE &t, GTVector<GTVector<GFTYPE>*> &u, const GFTYPE &dt)> 
        stcallback = steptop_callback; // set tmp function with proper signature for...
    eqn_impl->set_bdy_update_callback(fcallback);
    eqn_impl->set_steptop_callback(stcallback);
    eqn_impl->set_nu(nu_);

    // Create the mixer (to update forcing)
    EH_MESSAGE("main: create mixer...");
    MixBasePtr pMixer;
    create_mixer(ptree, pMixer);

    // Initialize state:
    EH_MESSAGE("main: Initializing state...");

    for ( GSIZET j=0; j<u_.size(); j++ ) {
      sprintf(stmp, "u%da", j+1);
      savars.push_back(stmp);
    }
    for ( GSIZET j=0; j<c_.size(); j++ ) {
      sprintf(stmp, "c%d", j+1);
      scvars.push_back(stmp);
      sprintf(stmp, "diff%d", j+1);
      sdvars.push_back(stmp);
    }
    if ( itindex == 0 ) { // start new run
      icycle = 0;
      t      = 0.0; 
      compute_analytic(*grid_, t, ptree, u_);
      compute_analytic(*grid_, t, ptree, ua_); // for setting error norm
    }
    else {              // restart run
      gio_restart(ptree, 0, u_, p, icycle, t, comm_);
    }

    // Create the observers: 
    EH_MESSAGE("main: create observers...");
    std::shared_ptr<std::vector<std::shared_ptr<ObserverBase<MyTypes>>>> pObservers(new std::vector<std::shared_ptr<ObserverBase<MyTypes>>>());
    create_observers(ptree, icycle, t, eqn_base, pObservers);
    for ( GSIZET j=0; j<pObservers->size(); j++ ) (*pObservers)[j]->set_tmp(utmp_);

    // Create integrator:
    EH_MESSAGE("main: create integrator...");
    auto pIntegrator = IntegratorFactory<MyTypes>::build(tintptree, eqn_base, pMixer, pObservers, *grid_);
    pIntegrator->get_traits().cycle = icycle;


    // Do time integration (output included
    // via observer(s)):
    EH_MESSAGE("main: do time stepping...");
    GPP(comm_,"main: do time stepping...");
    GTimerReset();
    GTimerStart("time_loop");

    pIntegrator->time_integrate(t, uf_, ub_, u_);

    GTimerStop("time_loop");
    GPP(comm_,"main: time stepping done.");

#if 1
    // Compute analytic solution, do comparisons:
    EH_MESSAGE("main: gather errors...");
    GPP(comm_,"main: gather errors...");
    GFTYPE           dxmin, lmin;
    GTVector<GFTYPE> lnorm(3), gnorm(3), maxerror(3);
    GTVector<GFTYPE> nnorm(nsolve_);
    GString sdir = ".";
    
    maxerror = 0.0;
    lnorm    = 0.0;  
    nnorm    = 1.0;


    cout << "main: t_end=" << t << endl;
    tt = 0.0;
    compute_analytic(*grid_, tt, ptree, ua_); // analyt soln at t=0
    for ( GSIZET j=0; j<nsolve_; j++ ) { // local errors
     *utmp_[1] = *ua_[j]; utmp_[1]->pow(2);
      nnorm[j] = grid_->integrate(*utmp_  [1],*utmp_[0]) ; // L2 norm of analyt soln at t=0
      nnorm[j] = nnorm[j] > std::numeric_limits<GFTYPE>::epsilon() ? nnorm[j] : 1.0;
      cout << "main: nnorm[" << j << "]=" << nnorm[j] << endl;
    }
    
    GTVector<GINT> istate(nsolve_);
    GTVector<GINT> cstate(c_.size());
    for ( GSIZET j=0; j<nsolve_; j++ ) istate[j] = j;
    for ( GSIZET j=0; j<c_.size(); j++ ) cstate[j] = j;

#if 1
    compute_analytic(*grid_, t, ptree, ua_); // analyt soln at t
    GIOTraits iot;
    iot.nelems = grid_->nelems();
    iot.gtype  = grid_->gtype();
    iot.porder.resize(1,GDIM);
    for ( GSIZET j=0; j<GDIM; j++ ) iot.porder(0,j) = gbasis[j]->getOrder();
    gio_write_state(iot, *grid_, ua_, istate, savars, comm_);
    for ( GSIZET j=0; j<GDIM; j++ ) 
    gio_write_state(iot, *grid_, c_, cstate, scvars, comm_);
    for ( GSIZET j=0; j<nsolve_; j++ ) { 
      *utmp_[j] = *u_[j] - *ua_[j];
    }
    gio_write_state(iot, *grid_, utmp_, istate, sdvars, comm_);
#endif

    for ( GSIZET j=0; j<nsolve_; j++ ) { //local errors
     *utmp_[0] = *u_[j] - *ua_[j];
#if 0
      GPP(comm_,"main: diff=u-ua[" << j << "]=" << *utmp_[0]);
#endif
     *utmp_[1]  = *utmp_[0]; utmp_[1]->abs();
     *utmp_[2]  = *utmp_[0]; utmp_[2]->pow(2);
      lnorm[0]  = utmp_[0]->infnorm (); // inf-norm
      gnorm[1]  = grid_->integrate(*utmp_[1],*utmp_[0])/sqrt(nnorm[j]) ; // L1-norm
      gnorm[2]  = grid_->integrate(*utmp_[2],*utmp_[0]) ; // L2-norm
      // Accumulate to find global errors for this field:
      GComm::Allreduce(lnorm.data()  , gnorm.data()  , 1, T2GCDatatype<GFTYPE>() , GC_OP_MAX, comm_);
      gnorm[2] =  sqrt(gnorm[2]/nnorm[j]);
      // now find max errors of each type for each field:
      for ( GSIZET i=0; i<3; i++ ) maxerror[i] = MAX(maxerror[i],fabs(gnorm[i]));
    }

    dxmin = grid_->minnodedist();
    lmin  = grid_->minlength();
    if ( myrank == 0 ) 
      cout << "main: maxerror = " << maxerror << endl;
   
    GTVector<GSIZET> lsz(2), gsz(2);
    lsz[0] = grid_->nelems();
    lsz[1] = grid_->ndof();
    GComm::Allreduce(lsz.data(), gsz.data(), 2, T2GCDatatype<GSIZET>() , GC_OP_SUM, comm_);

    // Print convergence data to file:
    std::ifstream itst;
    std::ofstream ios;

    if ( myrank == 0 ) {
      itst.open("burgers_err.txt");
      ios.open("burgers_err.txt",std::ios_base::app);
  
      // Write header, if required:
      if ( itst.peek() == std::ofstream::traits_type::eof() ) {
        ios << "#ntasks" << "  ";
        ios << "ncyc"    << "  ";
        ios << "var_dt"  << "  ";
        for ( GSIZET j=0; j<GDIM; j++ ) ios << "p" << j+1 << "  ";
        ios << "num_elems    dx_min   EL_min     inf_err     L1_err      L2_err" << std::endl;
      }
      itst.close();
  
      ios << ntasks << "  " ;
      ios << pIntegrator->get_numsteps()  << "  ";
      ios << solver_traits.variabledt << "  ";
      for ( GSIZET j=0; j<GDIM; j++ ) ios << pstd[j] << "  ";
      ios << gsz[0] << "  " << dxmin << "  " << lmin
          << "  " << maxerror[0] << "  " << maxerror[1] << "  " << maxerror[2]
          << std::endl;
      ios.close();
    }

#endif

    do_bench("benchmark.txt", pIntegrator->get_numsteps());
 
#if defined(_G_USE_GPTL)
//  GPTLpr(myrank);
    GPTLpr_file("timings.txt");
    GPTLpr_summary();
#endif
    GTimerFinal();

    EH_MESSAGE("main: do shutdown...");
//  GComm::TermComm();
    GlobalManager::shutdown();
    GlobalManager::finalize();
    if ( grid_ != NULLPTR ) delete grid_;
    for ( GSIZET j=0; j<GDIM; j++ ) delete gbasis[j];
    for ( GSIZET j=0; j<uold_.size(); j++ ) delete uold_[j];
    for ( GSIZET j=0; j<utmp_.size(); j++ ) delete utmp_[j];
    for ( GSIZET j=0; j<u_   .size(); j++ ) delete u_   [j];
    for ( GSIZET j=0; j<ua_  .size(); j++ ) delete ua_  [j];
    for ( GSIZET j=0; j<ub_  .size(); j++ ) delete ub_  [j];
    for ( GSIZET j=0; j<uf_  .size(); j++ ) delete uf_  [j];

    return(0);

} // end, main


//**********************************************************************************
//**********************************************************************************
// METHOD : update_dirichlet
// DESC   : update/set Dirichlet vectors, ub
// ARGS   : t    : time
//          u    : current state
//          ub   : bdy vectors (one for each state element)
// RETURNS: none.
//**********************************************************************************
void update_dirichlet(const GFTYPE &t, GTVector<GTVector<GFTYPE>*> &u, GTVector<GTVector<GFTYPE>*> &ub)
{

  GFTYPE tt = t;
  GTVector<GTVector<GSIZET>> *igbdy = &grid_->igbdy();
  
  // If bc is time dependent, update it here.
  // Note: grid_ ptree, and ua_ are global:
  if ( (*igbdy)[GBDY_DIRICHLET].size() > 0 ) {
    compute_analytic(*grid_, tt, ptree,  ua_);
  }

  // ...GBDY_DIRICHLET:
  // Set from State vector, ua_:
  for ( GSIZET k=0; k<nsolve_; k++ ) { 
    for ( GSIZET j=0; j<(*igbdy)[GBDY_DIRICHLET].size(); j++ ) {
      (*ub[k])[j] = (*ua_[k])[(*igbdy)[GBDY_DIRICHLET][j]];
    }
  }

} // end of method update_dirichlet


//**********************************************************************************
//**********************************************************************************
// METHOD: compute_dirplnwave_burgers
// DESC  : Compute solution to (nonlinear) Burgers with planar or
//         cicularized N-wave initial conditions for use with GBDY_DIRICHLET bdys
//         Must use box grid.
// ARGS  : grid    : GGrid object
//         t       : time
//         ptree   : main property tree
//         ua      : return solution
//**********************************************************************************
void compute_dirplnwave_burgers(GGrid &grid, GFTYPE &t, const PropertyTree& ptree,  GTVector<GTVector<GFTYPE>*> &ua)
{
  GBOOL            bplanar=TRUE; // planar or circularized
  GBOOL            brot   =FALSE;
  GSIZET           i, j, idir, nxy;
  GFTYPE           A, K2, Re, r2, t0, tdenom;
  GFTYPE           efact, sum, tfact, xfact;
  GTVector<GFTYPE> K(GDIM), xx(GDIM), si(GDIM), sig(GDIM);
  GTPoint<GFTYPE>  r0(3), P0(3), gL(3);
  std::vector<GFTYPE> kprop;

  PropertyTree nwaveptree = ptree.getPropertyTree("init_nwave");
  PropertyTree boxptree   = ptree.getPropertyTree("grid_box");

  GTVector<GTVector<GFTYPE>> *xnodes = &grid_->xNodes();

  assert(grid.gtype() == GE_REGULAR && "Invalid element types");

  GTVector<GString> bc(6);
  bc[0] = boxptree.getValue<GString>("bdy_x_0");
  bc[1] = boxptree.getValue<GString>("bdy_x_1");
  bc[2] = boxptree.getValue<GString>("bdy_y_0");
  bc[3] = boxptree.getValue<GString>("bdy_y_1");
  bc[4] = boxptree.getValue<GString>("bdy_z_0");
  bc[5] = boxptree.getValue<GString>("bdy_z_1");
//assert(bc.multiplicity("GBDY_DIRICHLET") >= 2*GDIM
//    && "Dirichlet boundaries must be set on all boundaries");

  nxy = (*xnodes)[0].size(); // same size for x, y, z

  // From Whitham's book, in 1d:
  // u(x,t) = (x/t) [ 1 + sqrt(t/t0) (e^Re - 1)^-1 exp(x^2/(4 nu t))i ]^-1
  // were Re is 'Reynolds' number: Re = A / 2nu; can think of
  // A ~ U L scaling.
  // Set some parameters:
  r0.x1  = nwaveptree.getValue<GFTYPE>("x0"); 
  r0.x2  = nwaveptree.getValue<GFTYPE>("y0"); 
  r0.x3  = nwaveptree.getValue<GFTYPE>("z0"); 
  A      = nwaveptree.getValue<GFTYPE>("ULparm",1.0);
  Re     = nwaveptree.getValue<GFTYPE>("Re",100.0);
  t0     = nwaveptree.getValue<GFTYPE>("t0",0.04);
  bplanar= nwaveptree.getValue<GBOOL>("planar",TRUE);
  kprop  = nwaveptree.getArray<GFTYPE>("prop_dir");
  K      = kprop;
  K     *= 1.0/K.Eucnorm();

  K2     = 0.0 ; for ( GSIZET i=0; i<GDIM; i++ ) K2 += K[i]*K[i];

  // If prop direction has more than one component != 0. then
  // front is rotated (but still planar):
  for ( i=0, brot=TRUE; i<GDIM; i++ ) brot = brot && K[i] != 0.0 ;
  for ( i=0, idir=0; i<GDIM; i++ ) if ( K[i] > 0 ) {idir=i; break;}

  if ( t == 0.0 ) t = K2 * t0;
  nu_[0] = A/(2.0*Re); // set nu from Re


  for ( j=0; j<nxy; j++ ) {
    for ( i=0; i<GDIM; i++ ) {
      xx[i] = (*xnodes)[i][j] - r0[i];
      (*ua[i])[j] = 0.0;
    }
    if ( bplanar ) { // compute k.r for planar wave
      for ( i=0, sum=0.0; i<GDIM; i++ ) { 
        sum += K[i]*xx[i];
        xx[i] = 0.0;
      }
      xx[0] = sum;
    }
    for ( i=0, r2=0.0; i<GDIM; i++ ) r2 += xx[i]*xx[i];  

    // u(x,t) = (x/t) [ 1 + sqrt(t/t0) (e^Re - 1)^-1 exp(x^2/(4 nu t)) ]^-1
    tdenom = 1.0/(4.0*nu_[0]*t);
    tfact  = bplanar ? sqrt(t/t0): t/t0;
    efact  = tfact * exp(r2*tdenom) / ( exp(Re) - 1.0 );
    xfact  = 1.0 /( t * (  1.0 + efact ) );
    for ( i=0; i<GDIM; i++ ) (*ua[i])[j] = xx[i]*xfact;
    // dU1max = 1.0 / ( t * (sqrt(t/A) + 1.0) );
    // aArea  = 4.0*nu_[0]*log( 1.0 + sqrt(A/t) );
  }
  
} // end, compute_dirplnwave_burgers


//**********************************************************************************
//**********************************************************************************
// METHOD: compute_dirgauss_lump
// DESC  : Compute solution to pure advection equation with 
//         GBDY_DIRICHLET bcs, a Gaussian 'lump'. Must use box grid.
// ARGS  : grid    : GGrid object
//         t       : time
//         ptree   : main property tree
//         ua      : return solution
//**********************************************************************************
void compute_dirgauss_lump(GGrid &grid, GFTYPE &t, const PropertyTree& ptree,  GTVector<GTVector<GFTYPE>*>  &ua)
{
  GString          serr = "compute_dirgauss_lump: ";
  GBOOL            bContin;
  GINT             j, n;
  GFTYPE           argxp; 
  GFTYPE           nxy, sig0, u0;
  GTVector<GFTYPE> xx(GDIM), si(GDIM), sig(GDIM), ufact(GDIM);
  GTPoint<GFTYPE>  r0(3), P0(3);

  PropertyTree heatptree = ptree.getPropertyTree("init_lump");
  PropertyTree boxptree = ptree.getPropertyTree("grid_box");
  PropertyTree advptree  = ptree.getPropertyTree("burgers_traits");
  GBOOL doheat   = advptree.getValue<GBOOL>("doheat");
  GBOOL bpureadv = advptree.getValue<GBOOL>("bpureadv");

  GTVector<GTVector<GFTYPE>> *xnodes = &grid_->xNodes();

  assert(grid.gtype() == GE_REGULAR && "Invalid element types");

  std::vector<GFTYPE> cs;
  if ( bpureadv ) {
    cs = heatptree.getArray<GFTYPE>("adv_vel");
  }

  // Check bdy conditioins:
  GTVector<GString> bc(6);
  bc[0] = boxptree.getValue<GString>("bdy_x_0");
  bc[1] = boxptree.getValue<GString>("bdy_x_1");
  bc[2] = boxptree.getValue<GString>("bdy_y_0");
  bc[3] = boxptree.getValue<GString>("bdy_y_1");
  bc[4] = boxptree.getValue<GString>("bdy_z_0");
  bc[5] = boxptree.getValue<GString>("bdy_z_1");
  assert(bc.multiplicity("GBDY_DIRICHLET") >= 2*GDIM 
      && "Dirichlet boundaries must be set on all boundaries");

  nxy = (*xnodes)[0].size(); // same size for x, y, z
  
  r0.x1 = heatptree.getValue<GFTYPE>("x0"); 
  r0.x2 = heatptree.getValue<GFTYPE>("y0"); 
  r0.x3 = heatptree.getValue<GFTYPE>("z0"); 
  sig0  = heatptree.getValue<GFTYPE>("sigma"); 
  u0    = heatptree.getValue<GFTYPE>("u0"); 

  // Set velocity here. May be a function of time.
  // These point to components of state u_:
  for ( j=0; j<GDIM; j++ ) *c_ [j] = 0.0;
 
  if ( bpureadv ) {
     for ( j=0; j<GDIM; j++ ) *c_[j] = cs[j];
  }

  // Prepare for case where sig is anisotropic (for later, maybe):
  for ( GSIZET k=0; k<GDIM; k++ ) {
    sig  [k] = sqrt(sig0*sig0 + 4.0*t*nu_[0]); // constant viscosity only
    si   [k] = 1.0/(sig[k]*sig[k]);
    ufact[k] = u0*pow(sig0/sig[k],GDIM);
  }

  // Ok, return to assumption of isotropic nu: 
  for ( GSIZET j=0; j<nxy; j++ ) {
    // Note: following c t is actually Integral_0^t c(t') dt', 
    //       so if c(t) changes, change this term accordingly:
    for ( GSIZET i=0; i<GDIM; i++ ) xx[i] = (*xnodes)[i][j] - r0[i] - (*c_[i])[j]*t;
    argxp = 0.0;
    for ( GSIZET i=0; i<GDIM; i++ ) argxp += -pow(xx[i],2.0)*si[i];
   (*ua[0])[j] = ufact[0]*exp(argxp);
  }

//GPP(comm_,serr << "ua=" << *ua[0] );
  
} // end, compute_dirgauss_lump


//**********************************************************************************
//**********************************************************************************
// METHOD: compute_icosgauss
// DESC  : Compute solution to pure advection equation with 
//         a Gaussian, but using, advection vel. from Williamson et al., 
//         JCP 102:211 (1992). May include multiple identical Gaussian
//         'hills'
// ARGS  : grid    : GGrid object
//         t       : time
//         ptree   : main property tree
//         ua      : return solution
//**********************************************************************************
void compute_icosgauss(GGrid &grid, GFTYPE &t, const PropertyTree& ptree,  GTVector<GTVector<GFTYPE>*>  &ua)
{
  GString          serr = "compute_icosgauss: ";
  GBOOL            bContin;
  GINT             j, k, n, nlumps;
  GSIZET           nxy;
  GFTYPE           alpha, argxp; 
  GFTYPE           lat, lon;
  GFTYPE           x, y, z, r, s;
  GFTYPE           rad, u0 ;
  GFTYPE           vtheta, vphi;
  GFTYPE           tiny = std::numeric_limits<GFTYPE>::epsilon();
  GTPoint<GFTYPE>           rt(3);
  GTVector<GFTYPE>          xx(3);
  GTVector<GFTYPE>          si(4), sig(4), ufact(4);
  GTVector<GFTYPE>          latp(4), lonp(4);
  std::vector<GFTYPE>       c0(4), sig0(4);
  std::vector<GFTYPE>       lat0(4), lon0(4); // up to 4 lumps

  PropertyTree lumpptree = ptree.getPropertyTree("init_icosgauss");
  PropertyTree icosptree = ptree.getPropertyTree("grid_icos");
  PropertyTree advptree  = ptree.getPropertyTree("burgers_traits");
  GBOOL doheat   = advptree.getValue<GBOOL>("doheat");
  GBOOL bpureadv = advptree.getValue<GBOOL>("bpureadv");

  GTVector<GTVector<GFTYPE>> *xnodes = &grid_->xNodes();

  assert(grid.gtype() == GE_2DEMBEDDED && "Invalid element types");

  std::vector<GFTYPE> Omega;

  nxy = (*xnodes)[0].size(); // same size for x, y, z
  
  lat0  = lumpptree.getArray<GFTYPE>("latitude0"); // lat for each lump
  lon0  = lumpptree.getArray<GFTYPE>("longitude0"); // lon for each lump
  sig0  = lumpptree.getArray<GFTYPE>("sigma"); // sig for each lump
  c0    = lumpptree.getArray<GFTYPE>("c0");  // initial concentrations for each lump
  rad   = icosptree.getValue<GFTYPE>("radius"); 
  u0    = lumpptree.getValue<GFTYPE>("u0"); 
  alpha = lumpptree.getValue<GFTYPE>("alpha",0.0); 
  nlumps= lumpptree.getValue<GINT>("nlumps",1); 

  alpha *= PI/180.0;

  // Convert initial locations from degrees to radians,
  // & compute initial positions of lumps in Cart coords:
  for ( GSIZET k=0; k<nlumps; k++ ) {
    lat0[k] *= PI/180.0;
    lon0[k] *= PI/180.0;
  }

  // Set velocity here. Taken to be solid body rotation,
  // u = Omega X r, where rotation rate vector, Omega
  // is computed by rotating u0 x (0, 0, 1) about x axis by
  // an amount alpha. These point to components of state u_:
  if ( bpureadv ) {
    for ( k=0; k<nxy; k++ ) {
      x   = (*xnodes)[0][k]; y = (*xnodes)[1][k]; z = (*xnodes)[2][k];
      r   = sqrt(x*x + y*y + z*z);
      // Compute lat & long:
      lat = asin(z/r);
      lon = atan2(y,x);
      // u_lat = u_theta = -u0 sin(lon) sin(alpha)
      // u_lon = u_phi   =  u0 (cos(theta) cos(alpha) + sin(theta)cos(lon)sin(alpha) )
      (*utmp_[0])[k]  = -u0*sin(lon)*sin(alpha);
      (*utmp_[1])[k]  =  u0*(cos(lat)*cos(alpha) + sin(lat)*cos(lon)*sin(alpha) );
    }
    GMTK::vsphere2cart(*grid_, utmp_, GVECTYPE_PHYS, c_);
//  GMTK::constrain2sphere(*grid_, c_);
  }

  *ua[0] = 0.0;
  for ( GSIZET k=0; k<nlumps; k++ ) {

    // Allow different sigma/concentration for each lump:
    sig  [k] = sqrt(sig0[k]*sig0[k] + 4.0*t*nu_[0]); // constant viscosity only
    si   [k] = 1.0/(sig[k]*sig[k]);
    ufact[k] = c0[k]*pow(sig0[k]/sig[k],GDIM);

    // Find where lat/lon endpoint would be at t, if alpha=0:
    latp[k]  = lat0[k];
    lonp[k]  = lon0[k] + (u0/rad) * t;

    // Find where Cart endpoint would be at t, if alpha=0:
    rt[0] = rad*cos(latp[k])*cos(lonp[k]); 
    rt[1] = rad*cos(latp[k])*sin(lonp[k]); 
    rt[2] = rad*sin(latp[k]);
  
    // Now, rotate rt about x-axis by alpha to
    // find lat/lon of final position of lump:
    xx[0] = rt[0]; xx[1] = rt[1]; xx[2] = rt[2];
    if ( t > 0 ) {
      xx[1] =  cos(alpha)*rt[1] + sin(alpha)*rt[2];
      xx[2] = -sin(alpha)*rt[1] + cos(alpha)*rt[2];
    }
    latp[k]  = asin(xx[2]/rad);
    lonp[k]  = atan2(xx[1],xx[0]);

    rt[0] = rad*cos(latp[k])*cos(lonp[k]); 
    rt[1] = rad*cos(latp[k])*sin(lonp[k]); 
    rt[2] = rad*sin(latp[k]);
  
    for ( GSIZET j=0; j<nxy; j++ ) {
      // Note: following c t is actually Integral_0^t c(t') dt', 
      //       so if c becomes a function of t, this muct change:
  
      x   = (*xnodes)[0][j];
      y   = (*xnodes)[1][j];
      z   = (*xnodes)[2][j];
#if 1
      r   = sqrt(x*x + y*y + z*z);
      lat = asin(z/r);
      lon = atan2(y,x);
      // Compute arclength from point to where center _should_ be:
      s     = r*acos( sin(latp[k])*sin(lat) + cos(latp[k])*cos(lat)*cos(lon-lonp[k]) );
     (*ua[0])[j] += ufact[k]*exp(-s*s*si[k]);
#else
     argxp = pow(x-rt[0],2) + pow(y-rt[1],2) + pow(z-rt[2],2);
     (*ua[0])[j] += ufact[k]*exp(-argxp*si[k]);
#endif
    } // end, grid-point loop

  } // end, lump loop

//GPP(comm_,serr << "ua=" << *ua[0] );
  
} // end, compute_icosgauss


//**********************************************************************************
//**********************************************************************************
// METHOD: compute_icosdefgauss
// DESC  : Compute solution to pure advection equation with 
//         dual Gaussian lumps, but using adv. field from
//         Nair & Lauritzen JCP 229:8868 (2010). NOTE:
//         'Analytic' solution to this problem is not available;
//         must integrate long enough for lumps to end up back
//         where they started, then compate with the t=0 solution.
// ARGS  : grid    : GGrid object
//         t       : time
//         ptree   : main property tree
//         ua      : return solution
//**********************************************************************************
void compute_icosdefgauss(GGrid &grid, GFTYPE &t, const PropertyTree& ptree,  GTVector<GTVector<GFTYPE>*>  &ua)
{
  GString          serr = "compute_icosgauss: ";
  GBOOL            bContin;
  GINT             nlumps=2;
  GSIZET           nxy;
  GFTYPE           alpha, argxp; 
  GFTYPE           lat0, lon0;
  GFTYPE           lat1, lon1;
  GFTYPE           lat, lon;
  GFTYPE           x, y, z, r;
  GFTYPE           rad, sig0, u0;
  GFTYPE           vtheta, vphi, Tper;
  GFTYPE           tiny = std::numeric_limits<GFTYPE>::epsilon();
  GTVector<GFTYPE> xx(3), si(3), sig(3), ufact(3);
  GTPoint<GFTYPE>  rt(3);
  GTVector<GTPoint<GFTYPE>>  r0(nlumps);

  PropertyTree lumpptree = ptree.getPropertyTree("init_icosdefgauss");
  PropertyTree icosptree = ptree.getPropertyTree("grid_icos");
  PropertyTree advptree  = ptree.getPropertyTree("burgers_traits");
  GBOOL doheat           = advptree.getValue<GBOOL>("doheat");
  GBOOL bpureadv         = advptree.getValue<GBOOL>("bpureadv");

  GTVector<GTVector<GFTYPE>> *xnodes = &grid_->xNodes();

  assert(grid.gtype() == GE_2DEMBEDDED && "Invalid element types");

  std::vector<GFTYPE> Omega;

  nxy = (*xnodes)[0].size(); // same size for x, y, z
  
  rad   = icosptree.getValue<GFTYPE>("radius"); 
  lat0  = lumpptree.getValue<GFTYPE>("latitude0", 0.0); 
  lon0  = lumpptree.getValue<GFTYPE>("longitude0",PI/3.0); 
  lat1  = lumpptree.getValue<GFTYPE>("latitude1", 0.0); 
  lon1  = lumpptree.getValue<GFTYPE>("longitude1",-PI/3.0); 
  sig0  = lumpptree.getValue<GFTYPE>("sigma", 0.15); 
  u0    = lumpptree.getValue<GFTYPE>("u0", 1.0); 
  Tper  = lumpptree.getValue<GFTYPE>("Tper", 1.0); 

  // Convert from degrees to radians (for later):
  lat0 *= PI/180.0;
  lon0 *= PI/180.0;
  lat1 *= PI/180.0;
  lon1 *= PI/180.0;

  // Compute initial position of lumps in Cart coords:
  r0[0].x1 = rad*cos(lat0)*cos(lon0);
  r0[0].x2 = rad*cos(lat0)*sin(lon0);
  r0[0].x3 = rad*sin(lat0);

  r0[1].x1 = rad*cos(lat1)*cos(lon1);
  r0[1].x2 = rad*cos(lat1)*sin(lon1);
  r0[1].x3 = rad*sin(lat1);

  // Set velocity here. Is a function of time,
  // and must be called every timestep.
  // These point to components of state u_:
  if ( bpureadv ) {
    for ( GSIZET k=0; k<nxy; k++ ) {
      x   = (*xnodes)[0][k]; y = (*xnodes)[1][k]; z = (*xnodes)[2][k];
      r   = sqrt(x*x + y*y + z*z);
      // Colmpute lat & long:
      lat = asin(z/r);
      lon = atan2(y,x);
      (*utmp_[0])[k]  = 0.5*u0*    sin(    lon)   *cos    (lat)*cos(PI*t/Tper);
      (*utmp_[1])[k]  =     u0*pow(sin(0.5*lon),2)*sin(2.0*lat)*cos(PI*t/Tper);
    }
    GMTK::vsphere2cart(*grid_, utmp_, GVECTYPE_PHYS, c_);
//  GMTK::constrain2sphere(*grid_, c_);
  }

  // Prepare for case where sig is anisotropic (for later, maybe):
  for ( GSIZET k=0; k<3; k++ ) {
    sig  [k] = sqrt(sig0*sig0 + 4.0*t*nu_[0]); // constant viscosity only
    si   [k] = 1.0/(sig[k]*sig[k]);
    ufact[k] = pow(sig0/sig[k],GDIM);
  }

  *ua[0] = 0.0;
  for ( GSIZET k=0; k<nlumps; k++ ) { // add each lump to solution
    for ( GSIZET j=0; j<nxy; j++ ) {
      rt.x1  = (*xnodes)[0][j];
      rt.x2  = (*xnodes)[1][j];
      rt.x3  = (*xnodes)[2][j];
      argxp  = 0.0;
      for ( GSIZET i=0; i<3; i++ ) argxp += pow(rt[i] - r0[k][i],2.0);
     (*ua[0])[j] += ufact[0]*exp(-argxp*si[0]);
    }
  }

//GPP(comm_,serr << "ua=" << *ua[0] );
  
} // end, compute_icosdefgauss


//**********************************************************************************
//**********************************************************************************
// METHOD: compute_icosbell
// DESC  : Compute solution to pure advection equation with 
//         a cosine bell, from Williamson et al., JCP 102:211 (1992)
// ARGS  : grid    : GGrid object
//         t       : time
//         ptree   : main property tree
//         ua      : return solution
//**********************************************************************************
void compute_icosbell(GGrid &grid, GFTYPE &t, const PropertyTree& ptree,  GTVector<GTVector<GFTYPE>*>  &ua)
{
  GString          serr = "compute_icosbell: ";
  GBOOL            bContin;
  GINT             j, k, n;
  GSIZET           nxy;
  GFTYPE           alpha, argxp; 
  GFTYPE           lat0, lon0, R, s;
  GFTYPE           lat, lon;
  GFTYPE           x, y, z, r;
  GFTYPE           xt, yt, zt;
  GFTYPE           rad, sig0, u0, u;
  GTVector<GFTYPE> xx(3), si(3), sig(3), ufact(3);
  GTPoint<GFTYPE>  r0(3);

  PropertyTree lumpptree = ptree.getPropertyTree("init_icosbell");
  PropertyTree icosptree = ptree.getPropertyTree("grid_icos");
  PropertyTree advptree  = ptree.getPropertyTree("burgers_traits");
  GBOOL doheat   = advptree.getValue<GBOOL>("doheat");
  GBOOL bpureadv = advptree.getValue<GBOOL>("bpureadv");

  GTVector<GTVector<GFTYPE>> *xnodes = &grid_->xNodes();

  assert(grid.gtype() == GE_2DEMBEDDED && "Invalid element types");

  std::vector<GFTYPE> Omega;

  nxy = (*xnodes)[0].size(); // same size for x, y, z
  
  rad   = icosptree.getValue<GFTYPE>("radius"); 
  lat0  = lumpptree.getValue<GFTYPE>("latitude0", 0.0); 
  lon0  = lumpptree.getValue<GFTYPE>("longitude0",3.0*PI/2.0); 
  sig0  = lumpptree.getValue<GFTYPE>("sigma", 0.15); 
  alpha = lumpptree.getValue<GFTYPE>("alpha",0.0); 
  u0    = lumpptree.getValue<GFTYPE>("u0", 1.0); 

  sig0  = 0.0;

  // Compute initial position of lump in Cart coords:
  r0.x1 = rad * cos(lat0*PI/180.0) * cos(lon0*PI/180.0);
  r0.x2 = rad * cos(lat0*PI/180.0) * sin(lon0*PI/180.0);
  r0.x3 = rad * sin(lat0*PI/180.0);

  // Set velocity here. May be a function of time.
  // These point to components of state u_:
  for ( j=0; j<3; j++ ) *c_ [j] = 0.0;
 
  if ( bpureadv ) {
    for ( k=0; k<nxy; k++ ) {
      x   = (*xnodes)[0][k]; y = (*xnodes)[1][k]; z = (*xnodes)[2][k];
      r   = sqrt(x*x + y*y + z*z);
      // Colmpute lat & long:
      lat = asin(z/r);
      lon = atan2(y,x);
      // u_theta = -u0 sin(lon) sin(alpha)
      // u_lamda =  u0 (cos(theta) cos(alpha) + sin(theta)cos(lon)sin(alpha) )
      (*utmp_[0])[k]  = -u0*sin(lon)*sin(alpha);
      (*utmp_[1])[k]  = u0*(cos(lat)*cos(alpha) + sin(lat)*cos(lon)*sin(alpha) );
    }
    GMTK::vsphere2cart(*grid_, utmp_, GVECTYPE_PHYS, c_);
  }


  // Prepare for case where sig is anisotropic (for later, maybe):
  for ( GSIZET k=0; k<3; k++ ) {
    sig  [k] = sqrt(sig0*sig0 + 4.0*t*nu_[0]); // constant viscosity only
    si   [k] = 1.0/(sig[k]*sig[k]);
    ufact[k] = u0*pow(sig0/sig[k],GDIM);
  }

  for ( GSIZET j=0; j<nxy; j++ ) {
    // Note: following c t is actually Integral_0^t c(t') dt', 
    //       so if c(t) above changes, change this term accordingly:
    x   = (*xnodes)[0][j];
    y   = (*xnodes)[1][j];
    z   = (*xnodes)[2][j];
    r   = sqrt(x*x + y*y + z*z);
    lat = asin(z/r);
    lon = atan2(y,x);

    // move origin:
    xt   = r0.x1 + (*c_[0])[j]*t; 
    yt   = r0.x2 + (*c_[1])[j]*t; 
    zt   = r0.x3 + (*c_[2])[j]*t;
    lat0 = asin(zt/r);
    lon0 = atan2(yt,xt);
    R   = r/3.0;
    s   = r*acos( sin(lat0)*sin(lat) + cos(lat0)*cos(lat)*cos(lon-lon0) );
   (*ua[0])[j] = s < R ? 0.5*u0*(1.0+cos(PI*s/R)) : 0.0;
   (*ua[0])[j]*= ufact[0];
  }

//GPP(comm_,serr << "ua=" << *ua[0] );
  
} // end, compute_icosbell


//**********************************************************************************
//**********************************************************************************
// METHOD: compute_pergauss_lump
// DESC  : Compute solution to pure advection equation with 
//         GBDY_PERIODIC bcs, a Gaussian 'lump'. Must use box grid.
// ARGS  : grid    : GGrid object
//         t       : time
//         ptree   : main property tree
//         ua      : return solution
//**********************************************************************************
void compute_pergauss_lump(GGrid &grid, GFTYPE &t, const PropertyTree& ptree,  GTVector<GTVector<GFTYPE>*>  &ua)
{
  GBOOL            bContin;
  GSIZET           i, j, k, n;
  GFTYPE           iargp, iargm ;
  GFTYPE           isum , irat , prod;
  GFTYPE           sumn , eps;
  GFTYPE           nxy, pint, sig0, u0;
  GTVector<GFTYPE> f(GDIM), xx(GDIM), si(GDIM), sig(GDIM);
  GTPoint<GFTYPE>  r0(3), P0(3), gL(3);

  PropertyTree heatptree = ptree.getPropertyTree("init_lump");
  PropertyTree boxptree = ptree.getPropertyTree("grid_box");
  PropertyTree advptree  = ptree.getPropertyTree("burgers_traits");
  GBOOL doheat   = advptree.getValue<GBOOL>("doheat");
  GBOOL bpureadv = advptree.getValue<GBOOL>("bpureadv");

  GTVector<GTVector<GFTYPE>> *xnodes = &grid_->xNodes();

  assert(grid.gtype() == GE_REGULAR && "Invalid element types");

  eps = 1.0e-4*std::numeric_limits<GFTYPE>::epsilon();

  // Get periodicity length, gL:
  std::vector<GFTYPE> xyz0 = boxptree.getArray<GFTYPE>("xyz0");
  std::vector<GFTYPE> dxyz = boxptree.getArray<GFTYPE>("delxyz");
  P0 = xyz0; r0 = dxyz; gL = r0;

  std::vector<GFTYPE> cs;
  if ( bpureadv ) {
    cs = heatptree.getArray<GFTYPE>("adv_vel");
  }

  GTVector<GString> bc(6);
  bc[0] = boxptree.getValue<GString>("bdy_x_0");
  bc[1] = boxptree.getValue<GString>("bdy_x_1");
  bc[2] = boxptree.getValue<GString>("bdy_y_0");
  bc[3] = boxptree.getValue<GString>("bdy_y_1");
  bc[4] = boxptree.getValue<GString>("bdy_z_0");
  bc[5] = boxptree.getValue<GString>("bdy_z_1");
  assert(bc.multiplicity("GBDY_PERIODIC") >= 2*GDIM
      && "Periodic boundaries must be set on all boundaries");

  nxy = (*xnodes)[0].size(); // same size for x, y, z
  
  r0.x1 = heatptree.getValue<GFTYPE>("x0"); 
  r0.x2 = heatptree.getValue<GFTYPE>("y0"); 
  r0.x3 = heatptree.getValue<GFTYPE>("z0"); 
  sig0  = heatptree.getValue<GFTYPE>("sigma"); 
  u0    = heatptree.getValue<GFTYPE>("u0"); 

  // Set adv velocity components. Note:
  // First state elem is the scalar solution, and
  // the remainder are the velocity components:

  for ( j=0; j<GDIM; j++ ) {
    sig[j] = sqrt(sig0*sig0 + 4.0*t*nu_[0]);
    si [j] = 1.0/(sig[j]*sig[j]);
   *c_ [j] = 0.0;
  }

  // Set velocity here. May be a function of time.
  // These point to components of state u_:
  if ( bpureadv ) for ( j=0; j<GDIM; j++ ) *c_[j] = cs[j];

  for ( n=0; n<nxy; n++ ) {

    prod = 1.0;
    for ( k=0; k<GDIM; k++ ) {
      // Note: following c t is actually Integral_0^t c(t') dt', 
      //       so if c(t) changes, change this term accordingly:
      f [k]  = modf((*c_[k])[j]*t/gL[k],&pint);
//    f [k]  = (*c_[k])[n]*t/gL[k];
      xx[k]  = (*xnodes)[k][n] - r0[k] - f[k]*gL[k];

      isum    = 0.0;
      i       = 0;
      irat    = 1.0;
      while ( irat > eps ) { // inner sum
        iargp   = pow((xx[k]+i*gL[k]),2)*si[k];
        iargm   = pow((xx[k]-i*gL[k]),2)*si[k];
        sumn    = i==0 ? exp(-iargp) : exp(-iargp) + exp(-iargm);
        isum   += sumn;
        irat    = sumn / isum ;
        i++;
      }
      prod *= isum;
    }
    (*ua[0])[n] = u0*pow(sig0,GDIM)/pow(sig[0],GDIM)*prod;

  } // end, loop over grid points

  
} // end, compute_pergauss_lump



//**********************************************************************************
//**********************************************************************************
// METHOD: steptop_callback
// DESC  : Top-of-time-step callback ('backdoor') function. 
//         This function might, e.g. update the linear advection 
//         vel. components as a function of time, since they
//         are not solved for in a PDE, but are, rather, prescribed.
// ARGS  : 
//         t  : time
//         u  : current state
//         dt : time step
//**********************************************************************************
void steptop_callback(const GFTYPE &t, GTVector<GTVector<GFTYPE>*>  &u, const GFTYPE &dt)
{
  
#if 0
  GFTYPE tt = t;
  compute_analytic(*grid_, tt, ptree, ua_);
#endif

} // end, method steptop_callback


//**********************************************************************************
//**********************************************************************************
// METHOD: compute_analytic
// DESC  : Compute analytic solutions based on property tree
// ARGS  : grid    : GGrid object
//         t       : time
//         ptree   : main property tree
//         ua      : return solution
//**********************************************************************************
void compute_analytic(GGrid &grid, GFTYPE &t, const PropertyTree& ptree, GTVector<GTVector<GFTYPE>*>  &ua)
{
  PropertyTree advptree  = ptree.getPropertyTree("burgers_traits");
  GBOOL doheat   = advptree.getValue<GBOOL>("doheat");
  GBOOL bpureadv = advptree.getValue<GBOOL>("bpureadv");

  GString      sblock = ptree.getValue<GString>("init_block"); // name of initialization block
  GString      sbcs   = ptree.getValue<GString>("bdy_conditions"); // name of initialization block
  PropertyTree blockptree = ptree.getPropertyTree(sblock); // sub-block of main ptree describing initialization type

  if ( doheat  ) {
    
    if ( sblock == "init_lump" ) {
      if ( sbcs == "PERIODIC" ) {
      compute_pergauss_lump(grid, t, ptree, ua);
      } else { // is DIRICHLET
      compute_dirgauss_lump(grid, t, ptree, ua);
      }
    }
    else if ( sblock == "init_icosbell" ) {
      compute_icosbell(grid, t, ptree, ua);
    }
    else if ( sblock == "init_icosgauss" ) {
      compute_icosgauss(grid, t, ptree, ua);
    }
    else if ( sblock == "init_icosdefgauss" ) {
      compute_icosdefgauss(grid, t, ptree, ua);
    }
    else {
      assert(FALSE && "Invalid heat equation initialization specified");
    }
    return;
  } // end, heat equation inititilization types



  // Inititialize for pure advection:
  if ( bpureadv ) {
    if ( sblock == "init_lump"  ) {
      if ( sbcs == "PERIODIC" ) {
      compute_pergauss_lump(grid, t, ptree, ua);
      } else { // is DIRICHLET
      compute_dirgauss_lump(grid, t, ptree, ua);
      }
    }
    else if ( sblock == "init_icosgauss" ) {
      compute_icosgauss(grid, t, ptree, ua);
    }
    else if ( sblock == "init_icosdefgauss" ) {
      compute_icosdefgauss(grid, t, ptree, ua);
    }
    else if ( sblock == "init_icosbell" ) {
      compute_icosbell(grid, t, ptree, ua);
    }
    else {
      assert(FALSE && "Invalid pure adv equation initialization specified");
    }
    return;
  }

  // Inititialize for nonlinear advection:
  if ( sblock == "init_nwave" ) {
    compute_dirplnwave_burgers(grid, t, ptree, ua);
  }
  else {
    assert(FALSE && "Invalid Burgers equation initialization specified");
  }

} // end, compute_analytic


//**********************************************************************************
//**********************************************************************************
// METHOD: create_mixer
// DESC  : Create mixer from main ptree
// ARGS  : grid    : GGrid object
//         ggfx    : gather/scatter op, GGFX
//**********************************************************************************
void create_mixer(PropertyTree &ptree, MixBasePtr &pMixer)
{
    PropertyTree     mixptree;    // observer props 
    MixBase::Traits traits;

    GString smixer = ptree.getValue<GString>("default_mixer","none");
    if ( "none" == smixer ) {
      pMixer= std::make_shared<NullMixer<MyTypes>>(traits, *grid_);
    }
    else {
      mixptree = ptree.getPropertyTree(smixer);
      pMixer = MixerFactory<MyTypes>::build(mixptree, *grid_);
    }

} // end method create_mixer


//**********************************************************************************
//**********************************************************************************
// METHOD: create_observers
// DESC  : Create observer list from main ptree
// ARGS  : grid      : GGrid object
//         icycle    : initial icycle
//         time      : initial time
//         equation  : shared ptr to EqnBasePtr type
//         pObservers: gather/scatter op, GGFX
//**********************************************************************************
void create_observers(PropertyTree &ptree, GSIZET icycle, GFTYPE time,
std::shared_ptr<EqnBase> &equation,
std::shared_ptr<std::vector<std::shared_ptr<ObserverBase<MyTypes>>>> &pObservers)
{
    GINT    ivers;
    GSIZET  rest_ocycle;       // restart output cycle
    GSIZET  deltac;            // cycle interval
    GFTYPE  ofact;             // output freq in terms of restart output
    GFTYPE  deltat;            // time interval
    PropertyTree obsptree;     // observer props 
    GString dstr = "none";
    GString ptype;
    GString ctype;

    if ( bench_ ) return;

    std::vector<GString> default_obslist; default_obslist.push_back(dstr);
    std::vector<GString> obslist = ptree.getArray<GString>("observer_list",default_obslist);
    dstr = "constant";
    ptype = ptree.getValue<GString>("exp_order_type",dstr);

    // Tie cadence_type to restart type:
    obsptree = ptree.getPropertyTree("posixio_observer");
    ctype    = obsptree.getValue<GString>("cadence_type");
   
    // If doing a restart, set observer output
    // cycles to value relative to the restart output cycle:
    rest_ocycle = ptree.getValue <GSIZET>("restart_index");
    if ( "constant" == ptype ) ivers = 0;
    if ( "variable" == ptype ) ivers = 1;
    for ( GSIZET j=0; j<obslist.size(); j++ ) {
      if ( "none" != obslist[j] ) {
        obsptree = ptree.getPropertyTree(obslist[j]);
        // Set output version based on exp_order_type:
        if ( "constant" == ptype 
         && "posixio_observer" == obslist[j]  ) obsptree.setValue<GINT>("misc",ivers);

        ofact       = obsptree.getValue<GDOUBLE>("interval_freq_fact",1.0);
        deltat      = obsptree.getValue<GDOUBLE>("time_interval",0.01);
        deltac      = obsptree.getValue<GDOUBLE>("cycle_interval",1);
        // Set current time and output cycle so that observer can initialize itself
        // These should be hidden from the config file:
        if ( "posixio_observer" == obslist[j]  ) ofact = 1.0;
        obsptree.setValue <GSIZET>("start_ocycle",MAX(0.0,rest_ocycle*ofact));
        obsptree.setValue <GFTYPE>("start_time"  ,time);
        obsptree.setValue <GFTYPE>("time_interval", MAX(0.0,deltat/ofact));
        obsptree.setValue <GSIZET>("cycle_interval",MAX(1.0,deltac/ofact));
        obsptree.setValue<GString>("cadence_type",ctype);

        pObservers->push_back(ObserverFactory<MyTypes>::build(obsptree, equation, *grid_));
      }
    }

} // end method create_observers


//**********************************************************************************
//**********************************************************************************
// METHOD: do_bench
// DESC  : Do benchmark from GPTL timers
// ARGS  : fname     : filename
//         ncyc      : number time cycles to average over
//**********************************************************************************
void do_bench(GString fname, GSIZET ncyc)
{
    if ( !bench_ ) return;

#if defined(_G_USE_GPTL)

    GINT   myrank   = GComm::WorldRank(comm_);
    GINT   ntasks   = GComm::WorldSize(comm_);
    GINT   nthreads = 0;
    GFTYPE ttotal;
    GFTYPE tggfx;
    GFTYPE texch;
    std::ifstream itst;
    std::ofstream ios;
    GTVector<GSIZET> lsz(2), gsz(2);

    // Get global no elements and dof:
    lsz[0] = grid_->nelems();
    lsz[1] = grid_->ndof();
    GComm::Allreduce(lsz.data(), gsz.data(), 2, T2GCDatatype<GSIZET>() , GC_OP_SUM, comm_);
    if ( myrank == 0 ) {
      itst.open(fname);
      ios.open(fname,std::ios_base::app);
  
      // Write header, if required:
      if ( itst.peek() == std::ofstream::traits_type::eof() ) {
        ios << "#nelems"  << "  ";
        ios << "ndof"     << "  ";
        ios << "ntasks"   << "  ";
        ios << "nthreads" << "  ";
        ios << "ttotal"   << "  ";
        ios << "tggfx"    << "  ";
        ios << "texch"           ;
        ios << endl;
      }
      itst.close();

      GPTLget_wallclock("time_loop"     , 0,  &ttotal); ttotal /= ncyc;
      GPTLget_wallclock("ggfx_doop"     , 0,  &tggfx ); tggfx  /= ncyc;
      GPTLget_wallclock("ggfx_doop_exch", 0,  &texch ); texch  /= ncyc;
  
      ios << gsz[0]          << "   " ;
      ios << gsz[1]          << "   " ;
      ios << ntasks          << "   " ;
      ios << nthreads        << "   ";
      ios << ttotal          << "   ";
      ios << tggfx           << "   ";
      ios << texch                   ;
      ios << endl;

      ios.close();
    }
#endif

    return;

} // end method do_bench
