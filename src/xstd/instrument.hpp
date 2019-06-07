/*
 * instrument.hpp
 *
 *  Created on: May 27, 2019
 *      Author: bflynt
 */

#ifndef INCLUDE_XSTD_INSTRUMENT_HPP_
#define INCLUDE_XSTD_INSTRUMENT_HPP_


#include "xstd/instrument/profiler.hpp"
#include "xstd/instrument/tracer.hpp"

#define INSTRUMENT(mssg) PROFILE(mssg) TRACER(mssg)

#endif /* INCLUDE_XSTD_INSTRUMENT_HPP_ */
