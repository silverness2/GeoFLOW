/*
 * output_icos.cpp
 *
 *  Created on: Jun 20, 2019
 *      Author: bryan.flynt
 */


#include "output_icos.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

void output_icos(const int nstep, const IcosGrid& grid, const IcosSoln& soln){
	const int num_decimals = 3;

	std::string filename;
	std::ofstream myfile;

	// Write Height
	filename = "d.h_it." + std::to_string(nstep);
	myfile.open(filename);
	myfile.precision(num_decimals);
	myfile << std::fixed; // std::scientific;
	for(std::size_t i = 1; i <= grid.nip; ++i){
		myfile << std::setw(12) << grid.lon_v(i);
		myfile << std::setw(12) << grid.lat_v(i);
		myfile << std::setw(25) << soln(i,1);
		myfile << "\n";
	}
	myfile.close();

	filename = "d.vx_it." + std::to_string(nstep);
	myfile.open(filename);
	myfile.precision(num_decimals);
	myfile << std::fixed; // std::scientific;
	for(std::size_t i = 1; i <= grid.nip; ++i){
		myfile << std::setw(12) << grid.lon_v(i);
		myfile << std::setw(12) << grid.lat_v(i);
		myfile << std::setw(25) << soln(i,2);
		myfile << "\n";
	}
	myfile.close();

	filename = "d.vy_it." + std::to_string(nstep);
	myfile.open(filename);
	myfile.precision(num_decimals);
	myfile << std::fixed; // std::scientific;
	for(std::size_t i = 1; i <= grid.nip; ++i){
		myfile << std::setw(12) << grid.lon_v(i);
		myfile << std::setw(12) << grid.lat_v(i);
		myfile << std::setw(25) << soln(i,3);
		myfile << "\n";
	}
	myfile.close();

	filename = "d.vz_it." + std::to_string(nstep);
	myfile.open(filename);
	myfile.precision(num_decimals);
	myfile << std::fixed; // std::scientific;
	for(std::size_t i = 1; i <= grid.nip; ++i){
		myfile << std::setw(12) << grid.lon_v(i);
		myfile << std::setw(12) << grid.lat_v(i);
		myfile << std::setw(25) << soln(i,4);
		myfile << "\n";
	}
	myfile.close();

}

