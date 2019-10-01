//==================================================================================
// Module       : gpdf1d.cpp
// Date         : 1/1/19 (DLR)
// Description  : GeoFLOW CGD analysis/postprocessing tool to compute
//                1d PDFs of stored binary output.
// Copyright    : Copyright 2019. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================

#include "gexec.h"
#include "gtypes.h"
#include <cstdio>
#include <unistd.h>
#include <iostream>
#include "gtvector.hpp"
#include "gtmatrix.hpp"
#include "gtmatrix.hpp"
#include "gtstat.hpp"
#include "gio.h"
#include "tbox/mpixx.hpp"
#include "tbox/global_manager.hpp"
#include "tbox/input_manager.hpp"


using namespace std;

GC_COMM      comm_ ;      // communicator


int main(int argc, char **argv)
{
    GBOOL   dolog=FALSE, bfixeddr=FALSE, bret;
    GINT    errcode=0, iopt; 
    GSIZET  nbins=500;
    GFTYPE  fmin=0.0, fmax=0.0;
    GString finput, foutput, spref="pdf1d";

    mpixx::environment env(argc,argv); // init GeoFLOW comm
    mpixx::communicator world;
    GlobalManager::initialize(argc,argv);
    GlobalManager::startup();

    comm_ = world; // need this for solver(s) & grid


    // Parse command line. ':' after char
    // option indicates that it takes an argument.
    while ((iopt = getopt(argc, argv, "b:g:l:p:u:h")) != -1) {
      switch (iopt) {
      case 'b': // handled by InputManager
          nbins = atoi(optarg);
          break;
      case 'g': // dolog?
          dolog = (GBOOL)atoi(optarg);
          break;
      case 'l': // lower dynamic range
          fmin     = atoi(optarg);
          bfixeddr = TRUE;
          break;
      case 'p': // set output file prefix
          spref = optarg;
          break;
      case 'u': // upper dynamic range
          fmax     = atoi(optarg);
          bfixeddr = TRUE;
          break;
      case ':': // missing option argument
          std::cout << argv[0] << ": option " << optopt << " requires an argument" << std::endl;
          exit(1);
          break;
      case '?':
      case 'h': // help
          std::cout << "usage: " << std::endl <<
          argv[0] << " [-h] [-b #bins] [-g dolog?] [-p file_prefix] [-l lower_dyn_range] [-u upper_dyn_range] file1 file2 ..." << std::endl;
          exit(1);
          break;
      default: // invalid option
          std::cout << argv[0] << ": option " << optopt << " invalid" << std::endl;
          exit(1);
          break;
      }
    }

    assert(argc > optind && "No files specified");


    GIOTraits                   giotraits;
    GTStat<GFTYPE>              gstat(nbins, comm_);
    GTVector<GFTYPE>            u; 

    // Currently we only know about GIO, so we assume
    // this is what the data is written in:
    giotraits.wfile  = 2048;
    giotraits.wtask  = 5;
    giotraits.wtime  = 6;
    giotraits.dir    = ".";

    // Process each file specified:
    for ( ; optind < argc; optind++ ) {
      finput  = argv[optind];
      foutput = spref + "_" + finput;

      // read in data
      gio_read(giotraits, finput, u);
      gstat.dopdf1d(u, bfixeddr, fmin, fmax, dolog, foutput); 
    }


    GlobalManager::shutdown();
    GlobalManager::finalize();


    exit(0);


} // end, main
