//==================================================================================
// Module       : gtest_hack.cpp
// Date         : 10/21/20 (DLR)
// Description  : GeoFLOW mini-app for tensor-producs performance 
//                improvement
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : 
//==================================================================================

#include "gexec.h"
#include "gtypes.h"
#include <cstdio>
#include <unistd.h>
#include <iostream>
#include "gptl.h"
#include <memory>
#include <cstdlib>
#include <cassert>
#include <random>
#include "gcomm.hpp"
#include "ggfx.hpp"
#include "gllbasis.hpp"
#include "gmass.hpp"
#include "gmorton_keygen.hpp"
#include "ggrid_factory.hpp"
#include "gmtk.hpp"
#include "tbox/property_tree.hpp"
#include "tbox/mpixx.hpp"
#include "tbox/global_manager.hpp"
#include "tbox/input_manager.hpp"
#include "gtools.h"
#include "gio.h"

using namespace geoflow::tbox;
using namespace std;

#define   _DO_REFDERIV
#undef   _DO_REFDERIVW
#if defined(_DO_REFDERIVW) && defined(_DO_REFDERIV)
  #error "Cannot define both _DO_REFDERIVW AND _DO_REFDERIV"
#endif

GGrid *grid_ = NULLPTR;

typedef  GTVector<GTVector<GFTYPE>*> State;

int main(int argc, char **argv)
{

    GINT    errcode ;
    GINT    nc=GDIM; // no. coords
    GFTYPE  eps=10.0*std::numeric_limits<GFTYPE>::epsilon();
    GFTYPE  dnorm, told=0.0, tnew=0.0;
    GString sgrid;// name of JSON grid object to use
    GTVector<GINT>
            pvec;
    std::vector<GINT> 
            pstd(GDIM);  
    GC_COMM comm = GC_COMM_WORLD;

    // Initialize comm:
    GComm::InitComm(&argc, &argv);
    GINT myrank  = GComm::WorldRank();
    GINT nprocs  = GComm::WorldSize();

    // Read main prop tree; may ovewrite with
    // certain command line args:
    PropertyTree ptree;     // main prop tree
    PropertyTree gridptree; // grid prop tree
    PropertyTree polyptree; // poly_test prop tree

    /////////////////////////////////////////////////////////////////
    /////////////////////// Initialize system ///////////////////////
    /////////////////////////////////////////////////////////////////
    ptree.load_file("hack_input.jsn");       // main param file structure

    // Create other prop trees for various objects:
    sgrid       = ptree.getValue<GString>("grid_type");
    pstd        = ptree.getArray<GINT>("exp_order");
    gridptree   = ptree.getPropertyTree(sgrid);
    polyptree   = ptree.getPropertyTree("poly_test");
    nc          = GDIM;
    assert(sgrid == "grid_box"); // use only Cartesian grids
    assert(pstd  >= GDIM      ); 

    pvec.resize(pstd.size(()); pvec = pstd; pvec.range(0,GDIM-1);
    
    
    // Set GTPL options:
#if defined(_G_USE_GPTL)
    // Set GTPL options:
    GPTLsetoption (GPTLcpu, 1);
    GPTLsetoption (GPTLsync_mpi, 1);
#endif
    // Initialize timer:
    GTimerInit();

    // Create basis:
    GTVector<GNBasis<GCTYPE,GFTYPE>*> gbasis(GDIM);
    for ( auto k=0; k<GDIM; k++ ) {
      gbasis [k] = new GLLBasis<GCTYPE,GFTYPE>(pstd[k]);
    }
    
    // Create grid:
    grid_ = GGridFactory::build(ptree, gbasis, comm);


    /////////////////////////////////////////////////////////////////
    ////////////////////// Allocate arrays //////////////////////////
    /////////////////////////////////////////////////////////////////
    GPTLstart("do_gather_op");
    // Initialize gather/scatter operator:
    GGFX<GFTYPE> ggfx;
    init_ggfx(ptree, *grid_, ggfx);
    GPTLstop("do_gather_op");

    // Create state and tmp space:
    GTVector<GTVector<GFTYPE>*> utmp (3);
    GTVector<GTVector<GFTYPE>*> dunew(nc);
    GTVector<GTVector<GFTYPE>*> duold(nc);
    GTVector<GTVector<GFTYPE>*> da   (nc);
    GTVector<GTVector<GFTYPE>*> du   (nc);
    GTVector<GFTYPE>            diff, u;
    
    for ( auto j=0; j<utmp .size(); j++ ) utmp [j] = new GTVector<GFTYPE>(grid_->size());
    for ( auto j=0; j<duold.size(); j++ ) duold[j] = new GTVector<GFTYPE>(grid_->size());
    for ( auto j=0; j<dunew.size(); j++ ) dunew[j] = new GTVector<GFTYPE>(grid_->size());
    for ( auto j=0; j<da   .size(); j++ ) da   [j] = new GTVector<GFTYPE>(grid_->size());
    u   .resize(grid_->size());
    diff.resize(grid_->size());
    du = NULLPTR;

    /////////////////////////////////////////////////////////////////
    ////////////////////// Compute solutions/////////////////////////
    /////////////////////////////////////////////////////////////////

    // Initialize u: set p, q, r exponents;
    //   u = x^p y^q z^r:
    std::vector<GFTYPE> 
            vpoly = polyptree.getArray<GFTYPE>("poly");
    
    GINT    ncyc  = polyptree.getValue<GINT>("ncycles",100);
    GFTYPE  p, q, r, x, y, z;
    GTVector<GTVector<GFTYPE>> 
           *xnodes = &grid_->xNodes();   

    p = vpoly[0]; q = vpoly[1]; r = GDIM > 2 ? vpoly[2] : 0.0;

    // Set function, and analytic derivative:
    dnorm = 0.0;
    for ( auto j=0; j<(*xnodes)[0].size(); j++ ) {
      x = (*xnodes)[0][j];
      y = (*xnodes)[1][j];
      z = GDIM == 3 ? (*xnodes)[2][j] : 1.0;
        u     [j] = pow(x,p)*pow(y,q)*pow(z,r);
      (*da[0])[j] = p==0 ? 0.0 : p*pow(x,p-1)*pow(y,q)*pow(z,r);
      (*da[1])[j] = q==0 ? 0.0 : q*pow(x,p)*pow(y,q-1)*pow(z,r);
      if ( xnodes->size() > 2 ) (*da[2])[j] = r==0 ? 0.0 : r*pow(x,p)*pow(y,q)*pow(z,r-1);
    } // end, loop over grid points

    for ( auto j=0; j<nc; j++ ) dnorm = MAX(dnorm, da[j]->amax());

    // Compute numerical derivs of u in each direction, using
    // different methods:
    GTimerStart("old_deriv");
    for ( auto j=0; j<duold.size(); j++ ) {  
        for ( auto n=0; n<ncyc; n++ ) {
          grid_->deriv(u, j+1, *utmp[0], *duold[j]);
        }
    } 
    GTimerStop("new_deriv");

    GTimerStart("new_deriv");
    for ( auto j=0; j<dunew.size(); j++ ) {  
        for ( auto n=0; n<ncyc; n++ ) {
          grid_->deriv(u, j+1, *utmp[0], *dunew[j]);
        }
    } 
    GTimerStop("new_deriv");

    GTVector<GFTYPE> 
                      lnorm(2), gnorm(2), maxerrold(2), maxerrnew((2);
    GTVector<GFTYPE> *maxerr;
    std::string       smethod[2] = {"old", "new"};

#if defined(_G_USE_GPTL)
    GPTLget_wallclock("old_deriv"     , 0,  &told); told /= ncyc;
    GPTLget_wallclock("new_deriv"     , 0,  &tnew); tnew /= ncyc;
#endif

    /////////////////////////////////////////////////////////////////
    //////////////////////// Compute Errors /////////////////////////
    /////////////////////////////////////////////////////////////////
    maxerrold = 0.0; maxerrnew = 0.0;
    for ( auto n=0; n<1; n++ ) { // over old and new methods
      for ( auto i=0; i<duoldw.size(); i++ ) du[i] = duold[i]; 
      maxerr = &maxerrold;
      if  ( n == 1 ) {
        for ( auto i=0; i<dunew.size(); i++ ) du[i] = dunew[i];
        maxerr = &maxerrnew;
      }
      for ( auto j=0; j<dunew.size(); j++ ) { // errors in each direction
        diff     = (*da[j]) - (*du[n])[j];
       *utmp [0] = diff;                  // for inf-norm
       *utmp [1] = diff; utmp[1]->pow(2); // for L2 norm

        lnorm[0] = utmp[0]->infnorm (); 
        gnorm[1] = sqrt(grid_->integrate(*utmp[2],*utmp[3]))/dnorm;

        // Accumulate to find global inf-norm:
        GComm::Allreduce(lnorm.data()  , gnorm.data()  , 1, T2GCDatatype<GFTYPE>() , GC_OP_MAX, comm);
cout << "main: gnorm[" << j << "]=" << gnorm << endl;
        // Now find max errors of each type for each field:
        for ( auto i=0; i<maxerr->size(); i++ ) (*maxerr)[i] = MAX((*maxerr)[i],gnorm[i]);
        if ( (*maxerr)[1] > eps ) {
          std::cout << "main: -------------------------------------derivative FAILED : direction=" << j << " method: " << smethod[n]  << std::endl;
          errcode = 1;
        } else {
          std::cout << "main: -------------------------------------derivative OK : direction=" << j << " method: " << smethod[n] << std::endl;
          errcode = 0;
        }

      } // end, direction loop

    } // end, new-old method loop


    // Print convergence data to file:
    std::ifstream itst;
    std::ofstream ios;
    itst.open("deriv_err.txt");
    ios.open("deriv_err.txt",std::ios_base::app);

    // Write header, if required:
    if ( itst.peek() == std::ofstream::traits_type::eof() ) {
    ios << "# p      num_elems  ncyc  inf_err_old  L2_err_old  t_old   inf_err_new  L2_err_new   t_new" << std::endl;
    }
    itst.close();

    ios << pvec << 
        << "  " << grid_->ngelems() 
        <  "  " << ncyc
        << "  " << maxerrold[0] << "  " << maxerrold[1] << "  " << told
        << "  " << maxerrnew[0] << "  " << maxerrnew[1] << "  " << tnew
        << std::endl;
    ios.close();
 
#if defined(_G_USE_GPTL)
    GPTLpr_file("timings.txt");
//  GPTLpr(GComm::WorldRank(comm_));
//  GPTLpr(0);
    GPTLpr_summary();
#endif
    GTimerFinal();

    GComm::TermComm();

    return( errcode );

} // end, main

