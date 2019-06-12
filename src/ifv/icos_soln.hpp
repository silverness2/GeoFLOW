/*
 * icos_soln.hpp
 *
 *  Created on: Jun 11, 2019
 *      Author: bflynt
 */

#ifndef SRC_IFV_ICOS_SOLN_HPP_
#define SRC_IFV_ICOS_SOLN_HPP_


#include <array>
#include <cstdint>
#include <vector>


struct IcosSoln {

public:
	using Real    = double;

	void resize(const int N);

	std::vector<Real>               h;
	std::vector<std::array<Real,3>> velo;
	//std::vector<std::array<Real,6>> ASV_ana;
};




#endif /* SRC_IFV_ICOS_SOLN_HPP_ */
