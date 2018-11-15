/*
 * stepper_concept.hpp
 *
 *  Created on: Nov 14, 2018
 *      Author: bryan.flynt
 */

#ifndef SRC_GEOFLOW_ODEINT_STEPPER_CONCEPT_HPP_
#define SRC_GEOFLOW_ODEINT_STEPPER_CONCEPT_HPP_


namespace geoflow {
namespace odeint {


template<typename T>
class StepperConcept {

public:

	using StepperType = T;

	template<typename SysType, typename StateType, typename TimeType>
	void step(SysType& sys, StateType& u, const TimeType& t, const TimeType& dt){
		this->impl()->step(sys,u,t,dt);
	}

	template<typename SysType, typename StateType, typename TimeType>
	void step(SysType& sys, const StateType& uin, const TimeType& t, StateType& uout, const TimeType& dt){
		this->impl()->step(sys,uin,t,uout,dt);
	}

protected:
	StepperType* impl(){
		return static_cast<StepperType*>(this);
	}

	StepperType const* impl() const{
		return static_cast<StepperType const*>(this);
	}
};

} // namespace odeint
} // namespace geoflow


#endif /* SRC_GEOFLOW_ODEINT_STEPPER_CONCEPT_HPP_ */
