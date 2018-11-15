

add_library(odeint INTERFACE)

#target_sources(odeint
#	PRIVATE
#)

#target_compile_features(odeint
#	PUBLIC
#	PRIVATE
#)

#target_compile_options(odeint 
#	PRIVATE
#		"$<$<CXX_COMPILER_ID:AppleClang>:-Wno-comment>"
#		"$<$<BOOL:${OpenMP_FOUND}>:${OpenMP_CXX_FLAGS}>"
#		"$<$<BOOL:${MPI_CXX_FOUND}>:${MPI_CXX_FLAGS}>"	
#)

#target_include_directories(odeint
#	PUBLIC
#		"$<INSTALL_INTERFACE:${MY_INSTALL_INCDIR}>"    
#		"$<BUILD_INTERFACE:${MY_PROJECT_INCDIR}>"
#		"${Boost_INCLUDE_DIRS}"
#		"$<$<BOOL:${MPI_CXX_FOUND}>:${MPI_CXX_INCLUDE_PATH}>"
#)

#target_link_libraries(odeint 
#	PRIVATE
#		tbox
#		"${Boost_LIBRARIES}"
#		"$<$<BOOL:${OpenMP_FOUND}>:${OpenMP_CXX_LIBRARIES}>"
#		"$<$<BOOL:${MPI_CXX_FOUND}>:${MPI_CXX_LIBRARIES}>"
#)

# Install public headers
install(
	DIRECTORY "${MY_PROJECT_INCDIR}/geoflow"
	DESTINATION "${MY_INSTALL_INCDIR}/geoflow"  
	FILES_MATCHING PATTERN "*.hpp"
)

# Install Library
install(TARGETS odeint
		LIBRARY DESTINATION "${MY_INSTALL_LIBDIR}"
		ARCHIVE DESTINATION "${MY_INSTALL_LIBDIR}")
		
		
		
		
		