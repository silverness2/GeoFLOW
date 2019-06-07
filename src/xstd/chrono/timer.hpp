/*
 * timer.hpp
 *
 *  Created on: Apr 22, 2019
 *      Author: bryan.flynt
 */

#ifndef TIMER_HPP_
#define TIMER_HPP_

#include <cstdlib>
#include <ctime>
#include <chrono>

namespace xstd {


/**
 * @brief
 * Detect the highest resolution clock between the system and steady
 */
using maxres_sys_or_steady =
		std::conditional<
		std::chrono::system_clock::period::den <= std::chrono::steady_clock::period::den,
		std::chrono::system_clock, std::chrono::steady_clock>::type;

/**
 * @brief
 * Detect the highest resolution non-sleeping clock
 */
using maxres_nonsleep_clock =
		std::conditional<
		std::chrono::high_resolution_clock::is_steady,
		std::chrono::high_resolution_clock, maxres_sys_or_steady>::type;

/**
 * @brief Timer struct to time events
 *
 * Timer to provide an uniform timer for all events.
 */
class Timer final {
public:
	using clock_type = maxres_nonsleep_clock;
	using time_point = clock_type::time_point;
	using duration   = clock_type::duration;


	void start() {
		start_ = now();
	}

	void stop() {
		stop_ = now();
	}

	double elapsed_time() const {
		auto secs = std::chrono::duration<double, std::chrono::seconds>(stop_ - start_);
		return secs.count();
	}

	time_point now(){
		return clock_type::now();
	}

protected:
	time_point start_;
	time_point stop_;

	time_point now(){
		return clock_type::now();
	}
};



} /* namespace xstd */


#endif /* TIMER_HPP_ */
