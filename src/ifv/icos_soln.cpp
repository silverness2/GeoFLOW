/*
 * icos_soln.cpp
 *
 *  Created on: Jun 11, 2019
 *      Author: bflynt
 */

#include "icos_soln.hpp"

IcosSoln::IcosSoln(const Integer N){
	this->resize(N);
}

void
IcosSoln::resize(const Integer N){
	h.resize(N);
	velo.resize(N);
}




