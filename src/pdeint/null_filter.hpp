/*
 * null_filter.hpp
 *
 *  Created on: Sep. 14, 2020
 *      Author: d.rosenberg
 */

#ifndef SRC_PDEINT_NULL_FILTER_HPP_
#define SRC_PDEINT_NULL_FILTER_HPP_


#include "pdeint/filter_base.hpp"

namespace geoflow {
namespace pdeint {

/**
 * Do nothing filter implementation
 *
 * The NullFilter is the default filter when nothing else is
 * provided and does nothing when called.
 *
 * @see FilterBase
 */
template<typename FilterType>
class NullFilter : public FilterBase<FilterType> {

public:
	using Interface     = FilterBase<FilterType>;
	using State         = typename Interface::State;
	using StateComp     = typename Interface::StateComp;
	using Grid          = typename Interface::Grid;
	using Time          = typename Interface::Time;
        using Filter        = FilterType;
        using FilterBase    = FilterBase<Filter>;
        using FilterBasePtr = std::shared_ptr<FilterBase>;

        
        NullFilter(Grid& grid){
        }

protected:

	void apply_impl(const Time&, /* time */const StateComp& /* u */, const State& /* utmp */, StateComp& /* uout */){
		// Do nothing ...
	}

	void apply_impl(const Time&, /* time */const StateComp& /* u */, const State& /* utmp */){
		// Do nothing ...
	}

};


} // namespace pdeint
} // namespace geoflow



#endif /* SRC_PDEINT_NULL_FILTER_HPP_ */
