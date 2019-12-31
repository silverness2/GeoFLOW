/*
 * null_stirrer.hpp
 *
 *  Created on: Mar 27, 2019
 *      Author: bflynt, d.rosenberg
 */

#ifndef SRC_PDEINT_NULL_STIR_HPP_
#define SRC_PDEINT_NULL_STIR_HPP_


#include "pdeint/stirrer_base.hpp"

namespace geoflow {
namespace pdeint {

/**
 * Do nothing Stir implementation
 *
 * The NullMixer is the default stirrer when nothing else is
 * provided and does nothing when called.
 *
 * @see MixerBase
 */
template<typename EquationType>
class NullMixer: public MixerBase<EquationType> {

public:
	using Interface  = MixerBase<EquationType>;
	using State      = typename Interface::State;
	using Grid       = typename Interface::Grid;
	using Time       = typename Interface::Time;
        
        NullMixer(typename MixerBase<EquationType>::Traits &traits, Grid &grid){
        }

protected:

	void stir_impl(const Time& /* t */, const State& /* u */, State& /* uf */){
		// Do nothing ...
	}

};


} // namespace pdeint
} // namespace geoflow



#endif /* SRC_PDEINT_NULL_STIR_HPP_ */
