

# Create the "configure.hpp" file
configure_file(
	"configure.hpp.in"
    "configure.hpp" 
)

add_subdirectory(tbox)