#
# Finds the correct "inline" keyword for system
#


message("Search for Inline Command:")

message(STATUS "Try command = [inline]")
check_cxx_source_compiles("
typedef int foo_t;
static inline foo_t static_foo(){return 0;}
foo_t foo(){return 0;}
int main(int argc, char *argv[]){return 0;}
"
INLINE_WORKS)
if( ${INLINE_WORKS} )
	set(INLINE_KEYWORD_STRING inline)
	message(STATUS "Found inline command: ${INLINE_KEYWORD_STRING}")
	return()
endif ()


message(STATUS "Try command = [__inline]")
check_cxx_source_compiles("
typedef int foo_t;
static __inline foo_t static_foo(){return 0;}
foo_t foo(){return 0;}
int main(int argc, char *argv[]){return 0;}
"
INLINE_WORKS)
if( ${INLINE_WORKS} )
	set(INLINE_KEYWORD_STRING __inline)
	message(STATUS "Found inline command: ${INLINE_KEYWORD_STRING}")
	return()
endif ()

message(STATUS "Try command = [__inline__]")
check_cxx_source_compiles("
typedef int foo_t;
static __inline__ foo_t static_foo(){return 0;}
foo_t foo(){return 0;}
int main(int argc, char *argv[]){return 0;}
"
INLINE_WORKS)
if( ${INLINE_WORKS} )
	set(INLINE_KEYWORD_STRING __inline__)
	message(STATUS "Found inline command: ${INLINE_KEYWORD_STRING}")
	return()
endif ()


