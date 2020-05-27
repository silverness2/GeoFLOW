/*
 * dyn_A.hpp
 *
 *  Created on: Jun 11, 2019
 *      Author: bflynt
 *
 *  This code is a direct translation into C++
 *  from the Fortran code provided by Yonggang Yu.
 */

#ifndef SRC_IFV_DYN_A_HPP_
#define SRC_IFV_DYN_A_HPP_

#include "icos_grid.hpp"
#include "icos_soln.hpp"


void dyn_A(const int iswcase, const IcosGrid& grid, const IcosSoln& soln, IcosSoln& acc);


#endif /* SRC_IFV_DYN_A_HPP_ */
