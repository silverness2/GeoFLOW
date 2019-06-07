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



	// Entry Point to your application
	// -->

	// Call shutdown call backs
	GlobalManager::shutdown();

	// Finalize global (once per run)
	GlobalManager::finalize();

	return 0;
}
