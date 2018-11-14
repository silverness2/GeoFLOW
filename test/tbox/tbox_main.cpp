/*
 * tbox.cpp
 *
 *  Created on: Nov 13, 2018
 *      Author: bflynt
 */

#include "geoflow/tbox/pio.hpp"
#include "geoflow/tbox/mpixx.hpp"
#include "geoflow/tbox/global_manager.hpp"

using namespace geoflow::tbox;



int main(int argc, char* argv[]) {

	// Start Up MPI
	mpixx::environment env(argc,argv);
	mpixx::communicator world;

	GlobalManager::initialize(argc,argv);


	pio::perr << "Testing\n";


	GlobalManager::finalize();

	return 0;
}
