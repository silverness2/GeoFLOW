/*
 * integrator.ipp
 *
 *  Created on: Nov 27, 2018
 *      Author: bflynt
 */


namespace geoflow {
namespace pdeint {


template<typename ET>
Integrator<ET>::Integrator(const EqnBasePtr&   equation,
		                     const MixerBasePtr& mixer,
		                     const ObsBasePtr&   observer,
                                     Grid&               grid,
		                     const Traits&       traits) :
	cycle_(0), traits_(traits), eqn_ptr_(equation), mixer_ptr_(mixer), obs_ptr_(observer), grid_(&grid) {
	ASSERT(nullptr != eqn_ptr_);
	ASSERT(nullptr != mixer_ptr_);
	ASSERT(nullptr != obs_ptr_);
}

template<typename ET>
void Integrator<ET>::time_integrate( Time&       t,
                                               State&      uf,
                                               State&      ub,
		                               State&      u ){
	ASSERT(nullptr != eqn_ptr_);
	ASSERT(nullptr != obs_ptr_);
	using std::min;

        Time t0 = t;
        Size n  = traits_.cycle_end - traits_.cycle ;
        if      ( traits_.integ_type == INTEG_CYCLE ) {
          steps( t0 
               , traits_.dt
               , n
               , uf
               , ub
               , u
               , t);
        } 
        else if ( traits_.integ_type == INTEG_TIME  ) {
          time ( t
               , traits_.time_end
               , traits_.dt
               , uf
               , ub
               , u);
          t = traits_.time_end;
        }

}


template<typename ET>
void Integrator<ET>::time( const Time& t0,
		                     const Time& t1,
		                     const Time& dt,
		                     State&      uf,
		                     State&      ub,
		                     State&      u ){
	ASSERT(nullptr != eqn_ptr_);
	ASSERT(nullptr != mixer_ptr_);
	ASSERT(nullptr != obs_ptr_);
	using std::min;

	Time t = t0;
	do {

		// Limit dt if requested
		Time dt_eff = dt;
		this->init_dt(t,u,dt_eff);
		dt_eff = min(dt_eff, t1 - t);
		if(0 >= dt_eff) {
			EH_ERROR("Effective Time Step is 0");
		}

		// Call observer on solution
                for ( auto j = 0 ; j < this->obs_ptr_->size(); j++ ) 
		  (*this->obs_ptr_)[j]->observe(t,u,uf);
 
		// Take Step
		this->eqn_ptr_->step(t, u, uf, ub, dt_eff);
		t += dt_eff;

                ++cycle_;

                // Call mixer to upate forcing:
		this->mixer_ptr_->mix(t,u, uf);


	} while( t != t1 );

	// Call observer on solution at final time:
        for ( auto j = 0 ; j < this->obs_ptr_->size(); j++ ) 
          (*this->obs_ptr_)[j]->observe(t,u,uf);
 


}

template<typename ET>
void Integrator<ET>::steps( const Time&  t0,
		                      const Time&  dt,
			              const Size&  n,
		                      State&       uf,
		                      State&       ub,
		                      State&       u,
			              Time&        t ){
	ASSERT(nullptr != eqn_ptr_);
	ASSERT(nullptr != mixer_ptr_);
	ASSERT(nullptr != obs_ptr_);


	t = t0;
	for(Size i = 0; i < n; ++i, ++cycle_){

		// Limit dt if requested
		Time dt_eff = dt;
		this->init_dt(t,u,dt_eff);
		if(0 >= dt_eff) {
			EH_ERROR("Effective Time Step is 0");
		}

		// Call observer on solution at new t
                for ( auto j = 0 ; j < this->obs_ptr_->size(); j++ ) 
		  (*this->obs_ptr_)[j]->observe(t,u,uf);

		// Take Step
		this->eqn_ptr_->step(t, u, uf, ub, dt_eff);
		t += dt_eff;

                // Call mixer to upate forcing:
		this->mixer_ptr_->mix(t,u, uf);

	}
        // Call observer on solution at final time:
        for ( auto j = 0 ; j < this->obs_ptr_->size(); j++ ) 
	  (*this->obs_ptr_)[j]->observe(t,u,uf);


}

template<typename ET>
void Integrator<ET>::list( const std::vector<Time>& tlist,
		                     State&                   uf,
		                     State&                   ub,
	                             State&                   u ){
	ASSERT(nullptr != eqn_ptr_);
	ASSERT(nullptr != obs_ptr_);

	for(Size i = 0; i < tlist.size()-1; ++i){

		// Calculate dt
		Time dt = tlist[i+1] - tlist[i];

		// Step till next time
		this->time(tlist[i],tlist[i+1],dt,uf,ub,u);
	}

}


template<typename ET>
void Integrator<ET>::init_dt(const Time& t, State& u, Time& dt) const{
	using std::max;
	using std::min;

	// Get pointer to equation
	auto eqn_ptr = this->eqn_ptr_; // ->getEquationPtr();
	ASSERT(nullptr != eqn_ptr);

	// Make sure dt is between absolute limits
	dt = min(dt, traits_.dt_max);
	dt = max(dt, traits_.dt_min);

	// Make sure dt is between CFL limits
	if( eqn_ptr->has_dt() ){
		Time dt_cfl_1 = std::numeric_limits<Time>::max()/(1+traits_.cfl_max);
//              eqn_ptr->dt(t,u,dt_cfl_1);
                eqn_ptr->dt(t,u,dt);
//              dt = min(dt, dt_cfl_1*traits_.cfl_max);
//              dt = max(dt, dt_cfl_1*traits_.cfl_min);
	}

}


} // namespace pdeint
} // namespace geoflow
