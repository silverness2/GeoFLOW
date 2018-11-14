#
# Check for C++11x smart pointers
#

check_cxx_source_compiles("

#include <memory>

int main(){
  std::unique_ptr<int> uip;
  std::shared_ptr<int> sip;
  std::weak_ptr<int>   wip;
  return 0;
}"
C11_SMARTPTR_WORKS) 
