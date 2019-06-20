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
	const auto sz = soln.size();
	const auto sdata = soln.data();
	const auto adata = acc.data();
	for(std::size_t i = 0; i < sz; ++i){
		sdata[i] += wtrk[istage-1] * adata[i];
	}
}


void
RK4th_test_Var_increment(const int istage, const double dt, const IcosSoln& soln_in, const IcosSoln& acc, IcosSoln& soln_out){
	//const double wtrk[4] = {dt/6, dt/3, dt/3, dt/6};
	const double stagecoef[4] = {0, dt/2, dt/2, dt};
	const auto sz = soln_in.size();
	const auto in_data  = soln_in.data();
	const auto out_data = soln_out.data();
	const auto adata    = acc.data();
	for(std::size_t i = 0; i < sz; ++i){
		out_data[i] = in_data[i] + stagecoef[istage-1] * adata[i];
	}
}
