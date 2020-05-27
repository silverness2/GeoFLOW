/*
 * output_icos.hpp
 *
 *  Created on: Jun 20, 2019
 *      Author: bryan.flynt
 *
 *  This code is a direct translation into C++
 *  from the Fortran code provided by Yonggang Yu.
 */

#ifndef SRC_IFV_OUTPUT_ICOS_HPP_
#define SRC_IFV_OUTPUT_ICOS_HPP_

#include "icos_grid.hpp"
#include "icos_soln.hpp"

void output_icos(const int nstep, const IcosGrid& grid, const IcosSoln& soln);


#endif /* SRC_IFV_OUTPUT_ICOS_HPP_ */
