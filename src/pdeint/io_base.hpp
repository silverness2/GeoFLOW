/*
 * io_base.hpp
 *
 *  Created on: Jan 21, 2020 
 *      Author: d.rosenberg
 */

#ifndef SRC_PDEINT_IOBASE_HPP_
#define SRC_PDEINT_IOBASE_HPP_


#include <memory>
#include <vector>

namespace geoflow {
namespace pdeint {

/**
 * IO base to carry out state IO 
 *
 */
template<typename TypePack>
class IOBase {

public:
        using Types        = TypePack;
	using State        = typename TypePack::State;
	using Grid         = typename TypePack::Grid;
	using Value        = typename TypePack::Value;
        using Time         = typename TypePack::Time;
	using Size         = typename TypePack::Size;
	using StateInfo    = typename TypePack::StateInfo; // May contain time, time index, var name etc
      
	IOBase() = delete;

	/**
	 * Constructor to initialize everything needed to do IO
	 *
	 */
	IOBase(){}
	IOBase(const IOBase& I) = default;
	~IOBase() = default;
	IOBase& operator=(const IOBase& I) = default;

	/**
	 * Write state at t
	 *
	 * @param[in,out] t  Initial time at start, and final time
	 * @param[in,out] u  Current state values
	 */
	void write_state( std::string  filename,
                          StateInfo&   info,
                          State&       u){
               return this->read_state_impl(filename, info, u);
             }
        /**
	 * Read state at t
	 *
	 * @param[in,out] t  Initial time at start, and final time
	 * @param[in,out] u  Current and final equation state values
	 */
	void read_state ( std::string  filename,
                          StateInfo&   info,
                          State&       u ){
               return this->write_state_impl(filename, info, u);
             }

protected:
        virtual void write_state_impl(std::string filename,
                                      StateInfo&  info,
                                      State&      u ) = 0;
        virtual void read_state_impl (std::string  filename,
                                      StateInfo&   info,
                                      State&       u ) = 0;
};


} // namespace pdeint
} // namespace geoflow


#include "pdeint/integrator.ipp"

#endif /* SRC_PDEINT_INTEGRATOR_HPP_ */
