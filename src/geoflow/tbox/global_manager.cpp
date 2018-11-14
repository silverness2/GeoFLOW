/*
 * global_manager.cpp
 *
 *  Created on: Nov 14, 2018
 *      Author: bflynt
 */

#include "geoflow/tbox/global_manager.hpp"

#include "geoflow/tbox/assert.hpp"
#include "geoflow/tbox/error_handler.hpp"
#include "geoflow/tbox/input_manager.hpp"
#include "geoflow/tbox/mpixx.hpp"
#include "geoflow/tbox/pio.hpp"

namespace geoflow {
namespace tbox {

bool GlobalManager::s_initialized = false;
bool GlobalManager::s_started     = false;
/*
 *************************************************************************
 *
 * Initialize the package.
 * This routine performs the following tasks:
 *
 * (1)
 * (2) Initialize the parallel I/O routines
 * (3) Set new handler so that an error message is printed if new fails.
 *
 *************************************************************************
 */
void GlobalManager::initialize(int argc, char* argv[]){
   ASSERT(!s_initialized);

	// Set ErrorHandler to use MPI
	EH::setAbortHandler(&mpixx_abort_handler);

	// Set new error handler
	std::set_new_handler(EH::badNew);

	// Start of PIO
	pio::initialize(mpixx::communicator().rank());
	//pio::logOnlyNodeZero("log");

	// Read Input File
	InputManager::initialize(argc,argv);

   //StartupShutdownManager::initialize();

   s_initialized = true;
}


void GlobalManager::startup(){
   ASSERT(s_initialized);
   ASSERT(!s_started);
   //StartupShutdownManager::startup();
   s_started = true;
}


void GlobalManager::shutdown(){
   ASSERT(s_initialized);
   ASSERT(s_started);
   //StartupShutdownManager::shutdown();
   s_started = false;
}


void GlobalManager::finalize(){
   ASSERT(s_initialized);
   //StartupShutdownManager::finalize();

   // Finalize PIO Buffers
   pio::finalize();

   s_initialized = false;
}

bool GlobalManager::isInitialized(){
   return s_initialized;
}

bool GlobalManager::isStarted(){
   return s_started;
}

} // namespace tbox
} // namespace geoflow

