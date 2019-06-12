/*
 * main.cpp
 *
 *  Created on: Nov 13, 2018
 *      Author: bflynt
 */

#include "tbox/pio.hpp"
#include "tbox/mpixx.hpp"
#include "tbox/global_manager.hpp"

#include "ifv/dyn_A.hpp"
#include "ifv/icos_grid.hpp"
#include "ifv/rk4.hpp"
#include "ifv/sw_test_init.hpp"
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

	// Create the IcosGrid & Soln
	IcosGrid grid;

	auto ptree = InputManager::getInputPropertyTree();
	auto nlevel   = ptree.getValue<int>("nlevel",4);
	auto stagger  = ptree.getValue<char>("stagger",'A');
	auto filename = ptree.getValue<std::string>("filename","sph_var_G04.dat");
	grid.load(stagger,nlevel,filename);

	/*

	// Run function "sw_test_init"
	IcosSoln soln_rk(grid.NPTS);
	auto iswcase = ptree.getValue<int>("iswcase",0);
	auto alpha   = ptree.getValue<double>("alpha",0);
	sw_test_init(iswcase,alpha,grid,soln_rk);

	// Save Initial Condition

	// Run SWM Dynamics
	auto itsbeg = ptree.getValue<int>("itsbeg",0);
	auto ForcastLength = ptree.getValue<int>("ForcastLength",100);
	double dt = 1456.0/std::pow(2.0,(nlevel-3));    // choice 1600/x, NICAM G4 dt=728  original comparable with NICAM
	int itsend = ForcastLength * 86400.0 / dt;
	IcosSoln afv(grid.NPTS);
	IcosSoln soln_0(grid.NPTS);
	IcosSoln soln_temp(grid.NPTS);
	for(int iter = itsbeg; iter < itsend; ++iter){
		// Zero afv but not needed
		soln_0.h    = soln_rk.h;
		soln_0.velo = soln_rk.velo;
		dyn_A(grid, soln_rk, afv);
		RK4th_true_Var_increment(0,dt,soln_rk,afv);

		for(int istage = 1; istage < 4; ++istage){
			RK4th_test_Var_increment(istage,dt,soln_0,afv,soln_temp);
			dyn_A(grid, soln_temp, afv);
			RK4th_true_Var_increment(istage,dt,soln_rk,afv);
		}
	}

	pio::pout << itsend << std::endl;

	*/


	// Entry Point to your application
	// -->

	// Call shutdown call backs
	GlobalManager::shutdown();

	// Finalize global (once per run)
	GlobalManager::finalize();

	return 0;
}
