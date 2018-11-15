/*
 * Clock.h
 *
 *  Created on: Mar 14, 2013
 *      Author: flyntbt
 */

#ifndef CLOCK_H_
#define CLOCK_H_

/**
 * Class Clock serves as a single point of access for system clock
 * information.
 *
 * The Clock class must be initialized before being used.
 *
 * The TimeStamp values are an always increasing value of 
 * seconds from some starting point value which is undefined.
 * Therefore, you must difference the values from multiple calls of 
 * TimeStamp() to get the time in seconds between the two calls.
 *
 * @code
 * time_diff = wallTimeStamp() - OldTimeStamp;  
 * @endcode
 *
 * User and System times only work if the back end library supports
 * the behavior otherwise zeros are returned in there place. 
 *
 *
 * - USE_C_TIME
 *     OS:
 *     - All
 *
 *     Requires;
 *     - ctime header from C
 *
 *     Resolution (seconds):
 *     - Wall clock   = 1
 *     - User clock   = 0.01
 *     - System clock = N/A
 *
 *
 * - USE_BOOST_TIME
 *     OS:
 *     - All
 *
 *     Requires;
 *     - boost::chrono
 *     - boost::system 
 *
 *     Resolution (seconds):
 *     - Wall clock   = 0.01
 *     - User clock   = 0.01
 *     - System clock = 0.01
 *
 *
 * - USE_C11_TIME
 *     OS:
 *     - All
 *
 *     Requires;
 *     - std::chrono (C++11)
 *
 *     Resolution (seconds):
 *     - Wall clock   = 0.01
 *     - User clock   = 0.01
 *     - System clock = N/A
 *
 *
 * - USE_POSIX_TIME
 *     OS:
 *     - Linux
 *
 *     Requires;
 *     - sys/times.h
 *     - sys/resource.h
 *     - sys/unistd.h
 *
 *     Resolution (seconds):
 *     - Wall clock   = 0.001
 *     - User clock   = 0.001
 *     - System clock = 0.001
 *
 *
 */

#include "geoflow/configure.hpp"

#if defined USE_C_TIME
#include <ctime>

#elif defined USE_C11_TIME
#include <chrono>

#elif defined USE_POSIX_TIME
#include <sys/times.h>
#include <sys/resource.h>
#include <sys/unistd.h>

#elif defined USE_BOOST_TIME
#include "boost/chrono.hpp"

#else
#define USE_C_TIME  // Default
#include <ctime>

#endif

// STL Includes
#include <string>



namespace geoflow {
namespace tbox {


struct Clock {
  
  /**
   * Initialize system clock.
   */
  static void initialize();
  
  /**
   * Timestamp for wall clock time
   */
  static double wallTimeStamp();
  
  /**
   * Timestamp for user clock time
   */
  static double userTimeStamp();
  
  /**
   * Timestamp for system clock time
   */
  static double systemTimeStamp();
  
  /**
   * Get a String representation of date/time
   */ 
  static std::string getDateStamp();
  
  
private:
  
#if defined USE_C_TIME
  static time_t  mStartWall;  // Measures Wall Time
  static clock_t mStartUser;  // Measures User Time

#elif defined USE_BOOST_TIME
  static boost::chrono::process_real_cpu_clock::time_point   mStartWall;
  static boost::chrono::process_user_cpu_clock::time_point   mStartUser;
  static boost::chrono::process_system_cpu_clock::time_point mStartSystem;

#elif defined USE_C11_TIME
  static std::chrono::steady_clock::time_point mStartWall;  // Measures Wall Time
  static clock_t                               mStartUser;  // Measures User Time

#elif defined USE_POSIX_TIME
  static clock_t       mStartWall;  // Measures Wall Time
  static struct rusage mStartUsage; // Measures User/Sys time
  static struct tms    mStartTMS;   // TMS type

#endif

};

} // namespace tbox
} // namespace geoflow


#endif /* CLOCK_H_ */
