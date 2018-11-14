#
# Tests for posix timer header and compile
#


check_cxx_source_compiles("

#include <sys/times.h>
#include <sys/resource.h>
#include <sys/unistd.h>

int main(){
  clock_t       mStartWall;  // Measures Wall Time
  struct rusage mStartUsage; // Measures User/Sys time
  struct tms    mStartTMS;   // TMS type

  mStartWall = times(&mStartTMS);
  getrusage(RUSAGE_SELF,&mStartUsage);
  
  double val = static_cast<double>(sysconf(_SC_CLK_TCK));

  return 0;    
}
" POSIX_TIME_WORKS )
