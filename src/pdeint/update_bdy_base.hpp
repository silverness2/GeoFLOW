/*
 * io_base.hpp
 *
 *  Created on: June 25, 2020 
 *      Author: d.rosenberg
 */

#ifndef SRC_PDEINT_UPDATEBDYBASE_HPP_
#define SRC_PDEINT_UPDATEBDYBASE_HPP_


#include <memory>
#include <vector>

namespace geoflow {
namespace pdeint {

/**
 * Base to handle specification and time updates of bdy information 
 * (or possibly state due to bdy) 
 *
 */
template<typename TypePack>
class UpdateBdyBase {

public:
        using Types        = TypePack;
	using State        = typename Types::State;
	using StateInfo    = typename Types::StateInfo; // May contain time, time index, var name etc
	using Grid         = typename Types::Grid;
	using Ftype        = typename Types::Value;
        using Time         = typename Types::Time;
	using IBdyVol      = typename Types::IBdyVol;
	using TBdyVol      = typename Types::TBdyVol;

      
	UpdateBdyBase() = default;
	UpdateBdyBase(const UpdateBdyBase& I) = default;
	~UpdateBdyBase() = default;
	UpdateBdyBase& operator=(const UpdateBdyBase& I) = default;

#if 0
	/**
	 * spec bdy conditions 
	 *
	 * @param[in]     sptree : property tree block for specification
	 * @param[in]     grid   : grid object
	 * @param[in]     bdyid  : bdy id
	 * @param[in]     ibdy   : indirection array into volmume indicating bdy indices
	 * @param[in,out] tbdy   : bdy types for each ibdy element
	 */
	bool spec   (PropertyTree &sptree,
                     Grid         &grid, 
                     int           bdyid, 
                     IBdyVol      &ibdy, 
                     TBdyvol      &tbdy){
                        return this->spec_impl(sptree, grid, bdyid, ibdy, tbdy);
                     }
#endif


	/**
	 * Update bdy conditions with state at t
	 *
	 * @param[in,out] grid   : initial time at start, and final time
	 * @param[in,out] stinfo : StateInfo object
	 * @param[in,out] time   : current time
	 * @param[in,out] utmp   : tmp arrays
	 * @param[in,out] u      : current state array
	 * @param[in,out] ub     : bdy arrays for each state component
	 */
	bool update (Grid      &grid, 
                     StateInfo &stinfo, 
                     Time      &time, 
                     State     &utmp, 
                     State     &u, 
                     State     &ub){
                        return this->update_impl(grid, stinfo, time, utmp, u, ub);
                     }

protected:

#if 0
	virtual bool spec   (PropertyTree &sptree,
                     Grid         &grid, 
                     int           bdyid, 
                     IBdyVol      &ibdy, 
                     TBdyvol      &tbdy) { return true;}
#endif
                     
	virtual bool update_impl (
                     Grid      &grid, 
                     StateInfo &stinfo, 
                     Time      &time, 
                     State     &utmp, 
                     State     &u, 
                     State     &ub) = 0;

};


} // namespace pdeint
} // namespace geoflow



#endif /* SRC_PDEINT_UPDATEBDYBASE_HPP_ */
