/*
 * display_header.cpp
 *
 *  Created on: Oct 3, 2019
 *      Author: bryan.flynt
 */


#include "configure.hpp"
#include "tbox/pio.hpp"

namespace geoflow {
namespace tbox {

void
display_header(const std::string& code_name){

    pio::pout << code_name << "\n";
    pio::pout << "Version:  "  << PROJECT_VERSION << "\n";
    pio::pout << "GIT Hash: " << GIT_LONG_HASH << "\n";
    pio::pout << "GIT Time: " << GIT_COMMIT_TIME << "\n";
    pio::pout << "Compile Time: " << __DATE__ << " " << __TIME__ << "\n";
    pio::pout << std::endl;
}


} // namespace tbox
} // namespace geoflow


