/*
 * system_base.hpp
 *
 *  Created on: Nov 14, 2018
 *      Author: bryan.flynt
 */

#ifndef SRC_GEOFLOW_ODEINT_SYSTEM_CONCEPT_HPP_
#define SRC_GEOFLOW_ODEINT_SYSTEM_CONCEPT_HPP_


namespace geoflow {
namespace odeint {


template<typename T>
class SystemConcept {

public:

	using SystemType = T;

	template<typename StateType, typename TimeType>
	void get_dt(const StateType& u, TimeType& dt){
		impl()->get_dt(u,dt);
	}

	template<typename StateType, typename DerivType, typename TimeType>
	void explicit_dudt(const StateType& u, DerivType& dudt, const TimeType& t = TimeType()){
		impl()->explicit_dudt(u,dudt,t);
	}

	template<typename StateType, typename DerivType, typename TimeType>
	void implicit_dudt(const StateType& u, DerivType& dudt, const TimeType& t = TimeType()){
		impl()->implicit_dudt(u,dudt,t);
	}

	template<typename StateType, typename JacType, typename TimeType>
	void implicit_dfdu(const StateType& u, JacType& dfdu, const TimeType& t = TimeType()){
		impl()->implicit_dfdu(u,dfdu,t);
	}

protected:
	SystemType* impl(){
		return static_cast<SystemType*>(this);
	}

	SystemType const* impl() const{
		return static_cast<SystemType const*>(this);
	}
};


} // namespace odeint
} // namespace geoflow



#endif /* SRC_GEOFLOW_ODEINT_SYSTEM_CONCEPT_HPP_ */
