#
# Tests for static_assert and compile
#

check_cxx_source_compiles("
#include <cassert>
int main(){
  assert(true);
  return 0;    
}
" C11_RUNTIME_ASSERT_WORKS )
