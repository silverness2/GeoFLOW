/*
 * null_observer.hpp
 *
 *  Created on: Nov 27, 2018
 *      Author: bflynt
 */

#ifndef SRC_PDEINT_NULL_OBSERVER_HPP_
#define SRC_PDEINT_NULL_OBSERVER_HPP_


#include "pdeint/observer_base.hpp"

namespace geoflow {
namespace pdeint {

/**
 * Do nothing Observer implementation
 *
 * The NullObserver is the default observer when nothing else is
 * provided and does nothing when called.
 *
 * @see ObserverBase
 */
template<typename EquationType>
class NullObserver : public ObserverBase<EquationType> {

public:
	using Interface  = ObserverBase<EquationType>;
	using State      = typename Interface::State;
	using Time       = typename Interface::Time;

protected:

	void observe_impl(const Time& /* t */, const State& /* u */){
		// Do nothing ...
	}

};


} // namespace pdeint
} // namespace geoflow



#endif /* SRC_PDEINT_NULL_OBSERVER_HPP_ */