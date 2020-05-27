/*
 * sw_test_init.hpp
 *
 *  Created on: Jun 7, 2019
 *      Author: bryan.flynt
 *
 *  NOTE:
 *	This code is a C++ representation of Fortran code
 *	provided by Yonggang Yu (CIRES/NOAA) to solve the
 *	2-D Shallow Water Equations on a sphere.
 */

#ifndef SRC_IFV_SW_TEST_INIT_HPP_
#define SRC_IFV_SW_TEST_INIT_HPP_

#include "icos_grid.hpp"
#include "icos_soln.hpp"

void sw_test_init(const int iswcase, const double alpha, IcosGrid& grid, IcosSoln& soln);



#endif /* SRC_IFV_SW_TEST_INIT_HPP_ */
