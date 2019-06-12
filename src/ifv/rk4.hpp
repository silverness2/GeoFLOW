/*
 * rk4.hpp
 *
 *  Created on: Jun 12, 2019
 *      Author: bryan.flynt
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
