/*
 * cpu_time.hpp
 *
 *  Created on: May 23, 2019
 *      Author: bflynt
 */

#ifndef CHRONO_CPU_TIME_HPP_
#define CHRONO_CPU_TIME_HPP_

#include <chrono>

#ifdef __linux__

#include <sys/times.h>
#include <unistd.h>
//#include <time.h>       // for clock_gettime

#elif __APPLE__

#include <sys/time.h>  //for gettimeofday and timeval
#include <sys/times.h> //for times
#include <unistd.h>

#endif



namespace xstd{
namespace chrono {


class process_real_cpu_clock {
public:
	typedef std::chrono::nanoseconds                   duration;
	typedef duration::rep                              rep;
	typedef duration::period                           period;
	typedef std::chrono::time_point<process_real_cpu_clock> time_point;
	static constexpr bool is_steady =                  true;

	static inline time_point now() noexcept;
};

class process_user_cpu_clock {
public:
	typedef std::chrono::nanoseconds                   duration;
	typedef duration::rep                              rep;
	typedef duration::period                           period;
	typedef std::chrono::time_point<process_user_cpu_clock> time_point;
	static constexpr bool is_steady =                  true;

	static inline time_point now() noexcept;
};

class process_system_cpu_clock {
public:
	typedef std::chrono::nanoseconds                     duration;
	typedef duration::rep                                rep;
	typedef duration::period                             period;
	typedef std::chrono::time_point<process_system_cpu_clock> time_point;
	static constexpr bool is_steady =                    true;

	static inline time_point now() noexcept;
};

namespace detail {
inline std::chrono::nanoseconds::rep tick_factor(){
	long factor = ::sysconf( _SC_CLK_TCK );
	factor = static_cast<long>(1000000000l) / factor;
    return factor;
}
}

process_real_cpu_clock::time_point process_real_cpu_clock::now() noexcept{
    tms tm;
    clock_t c = ::times( &tm );
    return time_point(std::chrono::nanoseconds(c*detail::tick_factor()));
}

process_user_cpu_clock::time_point process_user_cpu_clock::now() noexcept{
    tms tm;
    clock_t c = ::times( &tm );
    return time_point(std::chrono::nanoseconds((tm.tms_utime + tm.tms_cutime)*detail::tick_factor()));
}

process_system_cpu_clock::time_point process_system_cpu_clock::now() noexcept{
    tms tm;
    clock_t c = ::times( &tm );
    return time_point(std::chrono::nanoseconds((tm.tms_stime + tm.tms_cstime)*detail::tick_factor()));
}


} /* namespace chrono */
} /* namespace xstd */


/*
 * real = Wall clock time
 * user = Cumulative time spent by all the CPUs during the computation
 * sys  = Cumulative time spent by all the CPUs during system-related tasks such as memory, I/O, etc.
 * user + sys = Actual CPU time the process used.
 *
 * Function                                 Type         Resolution
 * ----------------------------------------------------------------------
 * std::chrono::system_clock                real         clock_gettime(CLOCK_REALTIME)
 * std::chrono::steady_clock                user+sys     clock_gettime(CLOCK_MONOTONIC)
 * time()                                   real         Seconds
 * times()                                  user, sys    10 MicroSeconds (100Hz)
 * clock()                                  user+sys     MicroSeconds
 * clock_gettime(CLOCK_REALTIME)            real         NanoSeconds
 * clock_gettime(CLOCK_MONOTONIC)           user+sys     NanoSeconds
 * clock_gettime(CLOCK_PROCESS_CPUTIME_ID)  user+sys     NanoSeconds Sum of Threads
 * clock_gettime(CLOCK_THREAD_CPUTIME_ID)   user+sys     NanoSeconds Single Thread
 * getrusage()                              user, sys    MicroSeconds
 *
 */





#endif /* CHRONO_CPU_TIME_HPP_ */
