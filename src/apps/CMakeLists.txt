#
# === Unit Tests
#

# 
# Build list of all tests to create
#
set(exec_files  
        geoflow.cpp
)

#
# Build list of compiler options for tests
#
set(compiler_options "$<$<CXX_COMPILER_ID:AppleClang>:-Wno-comment>"
		             "$<$<BOOL:${OpenMP_FOUND}>:${OpenMP_CXX_FLAGS}>"
		             "$<$<BOOL:${MPI_CXX_FOUND}>:${MPI_CXX_FLAGS}>" )
		       	
# Set include headers unique to main only
set(include_headers "")

# Set libraries the test need
# - Gets library include paths from here
set(linked_libraries tbox)	
				             
#
# Loop over each building & installing
#	             
FOREACH(test_src ${exec_files})
	string(REPLACE ".cpp" "" test_name ${test_src})
	add_executable(${test_name} ${test_src})
	target_compile_options(${test_name} PRIVATE ${compiler_options})
	target_include_directories(${test_name} PRIVATE ${include_headers})
	target_link_libraries(${test_name} ${linked_libraries})
	install(TARGETS ${test_name} DESTINATION "${MY_INSTALL_BINDIR}")
ENDFOREACH(test_src)

