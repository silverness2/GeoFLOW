/*
 * stop_watch.hpp
 *
 *  Created on: Jun 4, 2019
 *      Author: bryan.flynt
 */

#ifndef INCLUDE_XSTD_CHRONO_STOP_WATCH_HPP_
#define INCLUDE_XSTD_CHRONO_STOP_WATCH_HPP_

#include "xstd/chrono/timer.hpp"

#include <cassert>

class StopWatch final {
public:
	StopWatch() : active_(false){}
	StopWatch(const StopWatch& other) = default;
	StopWatch(StopWatch&& other) = default;
	~StopWatch() = default;
	StopWatch& operator=(const StopWatch& other) = default;
	StopWatch& operator=(StopWatch&& other) = default;


	void reset(){
		elapsed_ = duration::zero();
		active_  = false;
	}

	void start(){
		assert(not active_);
		active_ = true;
		start_  = clock_type::now();
	}

	void stop(){
		elapsed_ += (clock_type::now() - start_);
		assert(active_);
		active_ = false;
	}

	double elapsed_seconds() const{
		if(active_){
			auto now = clock_type::now();
			elapsed_ += (now - start_);
			start_ = now;
		}
		return std::chrono::duration<double, std::chrono::seconds>(elapsed_).count();
	}

private:
	using clock_type = xstd::maxres_nonsleep_clock;
	using time_point = clock_type::time_point;
	using duration   = clock_type::duration;

	time_point start_;
	duration   elapsed_;
	bool       active_;
};


#endif /* INCLUDE_XSTD_CHRONO_STOP_WATCH_HPP_ */
