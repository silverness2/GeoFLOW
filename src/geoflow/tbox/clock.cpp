/*
 * Clock.cpp
 *
 *  Created on: Mar 14, 2013
 *      Author: flyntbt
 */

#include "geoflow/tbox/clock.hpp"

// Always need this for Date Stamp
#include <ctime>
/*
 *************************************************************************
 *                                                                       *
 * Initialize Static Variables                                           *
 *                                                                       *
 *************************************************************************
 */
#if defined USE_C_TIME
time_t  Clock::mStartWall;  // Measures Wall Time
clock_t Clock::mStartUser;  // Measures User Time 

#elif defined USE_BOOST_TIME
boost::chrono::process_real_cpu_clock::time_point   Clock::mStartWall;
boost::chrono::process_user_cpu_clock::time_point   Clock::mStartUser;
boost::chrono::process_system_cpu_clock::time_point Clock::mStartSystem;

#elif defined USE_C11_TIME
std::chrono::steady_clock::time_point mStartWall;  // Measures Wall Time
clock_t Clock::mStartUser;                         // Measures User Time

#elif defined USE_POSIX_TIME
clock_t Clock::mStartWall;           // Measures Wall Time
struct rusage Clock::mStartUsage;    // Measures User/Sys time
struct tms Clock::mStartTMS;         // TMS type

#endif

/*
 *************************************************************************
 *                                                                       *
 * Initialize clock.                                                     *
 *                                                                       *
 *************************************************************************
 */
void Clock::initialize(){
#if defined USE_C_TIME
  mStartWall = time(NULL);
  mStartUser = clock();

#elif defined USE_BOOST_TIME
  mStartWall   = boost::chrono::process_real_cpu_clock::now();
  mStartUser   = boost::chrono::process_user_cpu_clock::now();
  mStartSystem = boost::chrono::process_system_cpu_clock::now();

#elif defined USE_C11_TIME
  mStartWall = std::chrono::steady_clock::now();
  mStartUser = clock();

#elif defined USE_POSIX_TIME
  mStartWall = times(&mStartTMS);
  getrusage(RUSAGE_SELF,&mStartUsage);

#endif
}

/*
 *************************************************************************
 *                                                                       *
 * Timestamp the provided structures with current system clock readings. *
 *                                                                       *
 *************************************************************************
 */
double Clock::wallTimeStamp(){
#if defined USE_C_TIME
  return difftime(time(NULL),mStartWall);

#elif defined USE_BOOST_TIME  
  boost::chrono::duration<double> dt = 
  boost::chrono::process_real_cpu_clock::now() - mStartWall;
  return dt.count();

#elif defined USE_C11_TIME
  std::chrono::duration<double> dt;
  dt = std::chrono::steady_clock::now() - mStartWall;
  return dt.count();

#elif defined USE_POSIX_TIME
  tms end_wall;
  clock_t end_clock = times(&end_wall);
  return (end_clock - mStartWall)/double(sysconf(_SC_CLK_TCK));

#endif
}
  

double Clock::userTimeStamp(){
#if defined USE_C_TIME
  return (clock() - mStartUser)/static_cast<double>(CLOCKS_PER_SEC);

#elif defined USE_BOOST_TIME  
  boost::chrono::duration<double> dt = 
  boost::chrono::process_user_cpu_clock::now() - mStartUser;
  return dt.count();

#elif defined USE_C11_TIME
  return (clock() - mStartUser)/static_cast<double>(CLOCKS_PER_SEC);

#elif defined USE_POSIX_TIME
  rusage end_time;
  getrusage(RUSAGE_SELF,&end_time);
  double sec = (end_time.ru_utime.tv_sec + (1.0e-6 * end_time.ru_utime.tv_usec))
    - (mStartUsage.ru_utime.tv_sec + (1.0e-6 * mStartUsage.ru_utime.tv_usec));
  return sec;

#endif
}

double Clock::systemTimeStamp(){
#if defined USE_C_TIME
  return 0;

#elif defined USE_BOOST_TIME  
  boost::chrono::duration<double> dt =
    boost::chrono::process_system_cpu_clock::now() - mStartSystem;
  return dt.count();

#elif defined USE_C11_TIME
  return 0;

#elif defined USE_POSIX_TIME
  rusage end_time;
  getrusage(RUSAGE_SELF,&end_time);
  double sec = (end_time.ru_stime.tv_sec - mStartUsage.ru_stime.tv_sec);
  sec += 0.000001*(end_time.ru_stime.tv_usec - mStartUsage.ru_stime.tv_usec);
  return sec;

#endif
}


/*
 *************************************************************************
 *                                                                       *
 * Get the Date Stamp in a std::string format                            *
 * Returns as "Day Month #Day hr:min:sec year"
 *                                                                       *
 *************************************************************************
 */
std::string Clock::getDateStamp(){
  time_t t = time(NULL);
  return ctime(&t);
}
