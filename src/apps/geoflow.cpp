/*
 * main.cpp
 *
 *  Created on: Nov 13, 2018
 *      Author: bflynt
 */

#include "tbox/pio.hpp"
#include "tbox/mpixx.hpp"
#include "tbox/global_manager.hpp"

#include "ifv/icos_grid.hpp"
#include "tbox/input_manager.hpp"

#include <cmath>

using namespace geoflow::tbox;

int main(int argc, char* argv[]) {

	// Start Up MPI
	mpixx::environment env(argc,argv);
	mpixx::communicator world;

	// Initialize global (once per run)
	GlobalManager::initialize(argc,argv);

	// Call startup call backs
	GlobalManager::startup();

	// Create the IcosGrid
	IcosGrid grid;

	auto ptree = InputManager::getInputPropertyTree();
	auto nlevel   = ptree.getValue<int>("nlevel",4);
	auto stagger  = ptree.getValue<char>("stagger",'A');
	auto filename = ptree.getValue<std::string>("filename","sph_var_G04.dat");
	grid.load(stagger,nlevel,filename);

	// Run function "sw_test_init"

	// Save Initial Condition

	// Run SWM Dynamics
	auto itsbeg = ptree.getValue<int>("itsbeg",0);
	auto ForcastLength = ptree.getValue<int>("ForcastLength",100);
	double dt = 1456.0/std::pow(2.0,(nlevel-3));    // choice 1600/x, NICAM G4 dt=728  original comparable with NICAM
	int itsend = ForcastLength * 86400.0 / dt;
	for(int iter = itsbeg; iter < itsend; ++iter){

	}



	// Entry Point to your application
	// -->

	// Call shutdown call backs
	GlobalManager::shutdown();

	// Finalize global (once per run)
	GlobalManager::finalize();

	return 0;
}
