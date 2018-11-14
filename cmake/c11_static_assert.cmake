#
# Tests for static_assert and compile
#

check_cxx_source_compiles("

#include <cassert>

int main(){
  static_assert(true,\"Testing for C++11x static_assert\");
  return 0;    
}
" C11_STATIC_ASSERT_WORKS )
