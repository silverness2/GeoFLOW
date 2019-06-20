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
#include "ifv/output_icos.hpp"
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


	// Run function "sw_test_init"
	IcosSoln soln_rk;
	soln_rk.resize(grid.NPTS,4).reorder(IcosSoln::f_order).rebase(1);
	auto iswcase = ptree.getValue<int>("iswcase",0);
	auto alpha   = ptree.getValue<double>("alpha",0);
	sw_test_init(iswcase,alpha,grid,soln_rk);

	// Save Initial Condition
	if( grid.stagger == 'A' ){
		output_icos(0,grid,soln_rk);
	}

	// Run SWM Dynamics
	auto itsbeg = ptree.getValue<int>("itsbeg",1);
	auto ForcastLength = ptree.getValue<int>("ForcastLength",100);
	double dt = 1456.0/std::pow(2.0,(nlevel-3));    // choice 1600/x, NICAM G4 dt=728  original comparable with NICAM
	int itsend = ForcastLength * 86400.0 / dt;

	IcosSoln afv, soln_0, soln_temp;
	afv.resize(grid.NPTS,4).reorder(IcosSoln::f_order).rebase(1);
	soln_0.resize(grid.NPTS,4).reorder(IcosSoln::f_order).rebase(1);
	soln_temp.resize(grid.NPTS,4).reorder(IcosSoln::f_order).rebase(1);

	for(int iter = itsbeg; iter <= itsend; ++iter){
		// Zero afv but not needed
		soln_0.copy(soln_rk); // TODO: this cause a reallocate
		dyn_A(iswcase, grid, soln_rk, afv);
		RK4th_true_Var_increment(1,dt,soln_rk,afv);

		for(int istage = 2; istage <= 4; ++istage){
			RK4th_test_Var_increment(istage,dt,soln_0,afv,soln_temp);
			dyn_A(iswcase, grid, soln_temp, afv);
			RK4th_true_Var_increment(istage,dt,soln_rk,afv);
		}

	     //------------------------------------------------------
	     //  output every 12 or 6 hours     !  post time = dt * iter
	     //------------------------------------------------------

		if( grid.stagger == 'A' ){
			auto nstep = int( 86400 / 2 / dt );
			if( iter%nstep == 0 ){
				auto j = int(iter/nstep);
				output_icos(j,grid,soln_rk);
			}
		}

	} // end iter loop

	// Call shutdown call backs
	GlobalManager::shutdown();

	// Finalize global (once per run)
	GlobalManager::finalize();

	return 0;
}
