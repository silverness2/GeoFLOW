

# Create the "configure.hpp" file
configure_file(
	"${MY_PROJECT_INCDIR}/geoflow/configure.hpp.in"
    "${MY_PROJECT_INCDIR}/geoflow/configure.hpp" 
)

add_subdirectory(tbox)