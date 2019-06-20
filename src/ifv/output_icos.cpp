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

void output_icos(const int nstep, const IcosGrid& grid, const IcosSoln& soln){
	std::string filename;
	std::ofstream myfile;

	// Write Height
	filename = "d.h_it." << nstep;
	myfile.open(filename);
	myfile.precision(2);
	myfile << std::fixed; // std::scientific;
	for(std::size_t i = 1; grid.nip; ++i){
		myfile << std::setw(12) << grid.lon_v(i);
		myfile << std::setw(12) << grid.lat_v(i);
		myfile << std::setw(25) << soln(i,1);
		myfile << "\n";
	}
	myfile.close();

	filename = "d.vx_it." << nstep;
	myfile.open(filename);
	myfile.precision(2);
	myfile << std::fixed; // std::scientific;
	for(std::size_t i = 1; grid.nip; ++i){
		myfile << std::setw(12) << grid.lon_v(i);
		myfile << std::setw(12) << grid.lat_v(i);
		myfile << std::setw(25) << soln(i,2);
		myfile << "\n";
	}
	myfile.close();

	filename = "d.vy_it." << nstep;
	myfile.open(filename);
	myfile.precision(2);
	myfile << std::fixed; // std::scientific;
	for(std::size_t i = 1; grid.nip; ++i){
		myfile << std::setw(12) << grid.lon_v(i);
		myfile << std::setw(12) << grid.lat_v(i);
		myfile << std::setw(25) << soln(i,3);
		myfile << "\n";
	}
	myfile.close();

	filename = "d.vz_it." << nstep;
	myfile.open(filename);
	myfile.precision(2);
	myfile << std::fixed; // std::scientific;
	for(std::size_t i = 1; grid.nip; ++i){
		myfile << std::setw(12) << grid.lon_v(i);
		myfile << std::setw(12) << grid.lat_v(i);
		myfile << std::setw(25) << soln(i,4);
		myfile << "\n";
	}
	myfile.close();

}

