#
# Check for C++17 <optional> header
#

check_cxx_source_compiles("
#include <cstdlib> // std::abort()
#include <optional>

std::optional<int> create(bool b){
	if( b )
		return 1;
	return {};
};

int main(){

	if( create( true ).value_or(0) != 1 ){
		std::abort();
	}
	
	if( create( false ).value_or(0) != 0 ){
		std::abort();
	}

  return 0;
}
"
C17_OPTIONAL_WORKS) 
