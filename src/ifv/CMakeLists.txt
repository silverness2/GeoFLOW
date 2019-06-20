

add_library(ifv STATIC "")

target_sources(ifv
	PRIVATE
	    dyn_A.cpp
		icos_grid.cpp
		icos_soln.cpp
		rk4.cpp
		sw_test_init.cpp
)

target_compile_features(ifv
	PUBLIC
		cxx_alias_templates
	PRIVATE
		cxx_nullptr	
)

target_compile_options(ifv 
	PRIVATE
		"$<$<CXX_COMPILER_ID:AppleClang>:-Wno-comment>"
		"$<$<BOOL:${OpenMP_FOUND}>:${OpenMP_CXX_FLAGS}>"
		"$<$<BOOL:${MPI_CXX_FOUND}>:${MPI_CXX_FLAGS}>"	
)

# Note: The configure.hpp file is within ${MY_INSTALL_SRCDIR}
# until it gets installed to ${MY_INSTALL_INCDIR} at the end.
target_include_directories(ifv
	PUBLIC  
		"$<INSTALL_INTERFACE:${MY_INSTALL_INCDIR}>"
		"$<BUILD_INTERFACE:${MY_PROJECT_INCDIR}>"
		"$<BUILD_INTERFACE:${MY_INSTALL_SRCDIR}>"
		"${Boost_INCLUDE_DIRS}"
		"$<$<BOOL:${MPI_CXX_FOUND}>:${MPI_CXX_INCLUDE_PATH}>"
)

target_link_libraries(ifv 
	PRIVATE
		"${Boost_LIBRARIES}"
		"$<$<BOOL:${OpenMP_FOUND}>:${OpenMP_CXX_LIBRARIES}>"
		"$<$<BOOL:${MPI_CXX_FOUND}>:${MPI_CXX_LIBRARIES}>"
)

# Install public headers
install(
	DIRECTORY "${MY_PROJECT_INCDIR}/ifv"
	DESTINATION "${MY_INSTALL_INCDIR}"  
	FILES_MATCHING PATTERN "*.hpp"
)

# Install Library
install(TARGETS ifv
		LIBRARY DESTINATION "${MY_INSTALL_LIBDIR}"
		ARCHIVE DESTINATION "${MY_INSTALL_LIBDIR}")
		
		
		
		
		