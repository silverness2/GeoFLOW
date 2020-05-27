/*
 * rk4.hpp
 *
 *  Created on: Jun 12, 2019
 *      Author: bryan.flynt

 *  NOTE:
 *	This code is a C++ representation of Fortran code
 *	provided by Yonggang Yu (CIRES/NOAA) to solve the
 *	2-D Shallow Water Equations on a sphere.
 */

#ifndef SRC_IFV_RK4_HPP_
#define SRC_IFV_RK4_HPP_

#include "icos_grid.hpp"
#include "icos_soln.hpp"


void
RK4th_true_Var_increment(const int istage, const double dt, IcosSoln& soln, const IcosSoln& acc);


void
RK4th_test_Var_increment(const int istage, const double dt, const IcosSoln& soln_in, const IcosSoln& acc, IcosSoln& soln_out);

#endif /* SRC_IFV_RK4_HPP_ */
