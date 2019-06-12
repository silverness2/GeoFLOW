/*
 * rk4.cpp
 *
 *  Created on: Jun 12, 2019
 *      Author: bryan.flynt
 */

#include "rk4.hpp"
#include "xstd/array.hpp"


void
RK4th_true_Var_increment(const int istage, const double dt, IcosSoln& soln, const IcosSoln& acc){
	const double wtrk[4] = {dt/6, dt/3, dt/3, dt/6};
	//const double stagecorf[4] = {0, dt/2, dt/2, dt};
	const auto sz = soln.h.size();
	for(std::size_t i = 0; i < sz; ++i){
		soln.h[i] += wtrk[istage] * acc.h[i];
	}
	for(std::size_t i = 0; i < sz; ++i){
		soln.velo[i] += wtrk[istage] * acc.velo[i];
	}
}


void
RK4th_test_Var_increment(const int istage, const double dt, const IcosSoln& soln_in, const IcosSoln& acc, IcosSoln& soln_out){
	//const double wtrk[4] = {dt/6, dt/3, dt/3, dt/6};
	const double stagecoef[4] = {0, dt/2, dt/2, dt};
	const auto sz = soln_in.h.size();
	for(std::size_t i = 0; i < sz; ++i){
		soln_out.h[i] = soln_in.h[i] + stagecoef[istage] * acc.h[i];
	}
	for(std::size_t i = 0; i < sz; ++i){
		soln_out.velo[i] = soln_in.velo[i] + stagecoef[istage] * acc.velo[i];
	}
}
