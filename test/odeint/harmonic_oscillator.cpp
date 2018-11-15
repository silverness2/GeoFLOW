/*
 * harmonic_oscillator.cpp
 *
 *  Created on: Nov 14, 2018
 *      Author: bryan.flynt
 */


#include "geoflow/odeint/system_concept.hpp"
#include "geoflow/odeint/stepper_concept.hpp"

#include <memory>

//
//
//
using namespace geoflow::odeint;


class HarmonicOscillator :
		public SystemConcept<HarmonicOscillator> {

private:
	static constexpr double gam = 0.15;

public:

	template<typename StateType, typename TimeType>
	void get_dt(const StateType& u, TimeType& dt){
		dt = 1.0;
	}

	template<typename StateType, typename DerivType, typename TimeType>
	void explicit_dudt(const StateType& u, DerivType& dudt, const TimeType& t){
		dudt[0] = +u[1];
		dudt[1] = -u[0] - gam*u[1];
	}

};




class BackwardEuler : StepperConcept<BackwardEuler>{

private:

public:

	template<typename SysType, typename StateType, typename TimeType>
	void step(SysType& sys, StateType& u, const TimeType& t, const TimeType& dt){
		StateType dudt(2); // How do we get and resize this ???
		sys.explicit_dudt(u,dudt,t);
		u[0] += dt * dudt[0];
		u[1] += dt * dudt[1];
	}

	template<typename SysType, typename StateType, typename TimeType>
	void step(SysType& sys, const StateType& uin, const TimeType& t, StateType& uout, const TimeType& dt){
		StateType dudt(2); // How do we get and resize this ???
		sys.explicit_dudt(uin,dudt,t);
		uout = uin;
		uout[0] += dt * dudt[0];
		uout[1] += dt * dudt[1];
	}

};


#include <vector>

int main() {

	HarmonicOscillator sys;
	BackwardEuler      stepper;


	const int N = 10;
	double t  = 0;
	double dt;
	std::vector<double> u(2, 0.0);

	// Take step
	for(int i = 0; i < N; ++i){
		sys.get_dt(u,dt);
		stepper.step(sys,u,t,dt);
		t += dt;
	}

	return 0;
}
