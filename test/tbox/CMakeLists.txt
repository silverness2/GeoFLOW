#
# === Unit Tests
#

# 
# Build list of all tests to create
#
set(test_source_files  mpixx.cpp
                       property_tree.cpp
                       tbox_main.cpp
)

#
# Build list of compiler options for tests
#
set(compiler_options "$<$<CXX_COMPILER_ID:AppleClang>:-Wno-comment>"
		             "$<$<BOOL:${OpenMP_FOUND}>:${OpenMP_CXX_FLAGS}>"
		             "$<$<BOOL:${MPI_CXX_FOUND}>:${MPI_CXX_FLAGS}>" )
		       
set(include_headers "${CMAKE_SOURCE_DIR}/src")
set(linked_libraries tbox)
		             
#
# Loop over each test building and adding to ctest
# - Build test name
#	             
FOREACH(test_src ${test_source_files})
	string(REPLACE ".cpp" "" test_name ${test_src})
	add_executable(${test_name} ${test_src})
	target_compile_options(${test_name} PRIVATE ${compiler_options})
	target_include_directories(${test_name} PRIVATE ${include_headers})
	target_link_libraries(${test_name} ${linked_libraries})
	add_test(${test_name}  ${test_name})
ENDFOREACH(test_src)


# Install test data
install(
	DIRECTORY "${MY_PROJECT_TESTDIR}/tbox/data"
	DESTINATION "${MY_INSTALL_TESTDIR}/tbox"  
	FILES_MATCHING PATTERN "*"
)
