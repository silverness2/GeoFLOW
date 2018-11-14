#
# Tests for chrono header and compile
#

check_cxx_source_compiles("

#include <chrono>

int main(){
  std::chrono::steady_clock::time_point mStartWall;
  mStartWall = std::chrono::steady_clock::now();
  return 0;    
}
" C11_CHRONO_WORKS )
