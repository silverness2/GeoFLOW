/*
 * icos_soln.hpp
 *
 *  Created on: Jun 11, 2019
 *      Author: bflynt
 */

#ifndef SRC_IFV_ICOS_SOLN_HPP_
#define SRC_IFV_ICOS_SOLN_HPP_


#include <cstdint>

#include "tbox/multi_array.hpp"

using IcosSoln = tbox::multi_array<double,2>;

/*
struct IcosSoln {

public:
	using Real    = double;
	using Integer = std::int32_t;
	using real_array_1 = tbox::multi_array<Real,1>;
	using real_array_2 = tbox::multi_array<Real,2>;
	static constexpr auto real_f_order_1 = real_array_1::f_order;
	static constexpr auto real_f_order_2 = real_array_2::f_order;
	static constexpr Integer ndim  = 3;

	IcosSoln(const Integer N);

	void resize(const Integer N);

	real_array_1    h;
	real_array_2 velo;
};
*/




#endif /* SRC_IFV_ICOS_SOLN_HPP_ */
