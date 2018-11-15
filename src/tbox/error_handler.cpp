

#include "tbox/error_handler.hpp"

#include <cstdlib> // exit()

#include "tbox/pio.hpp"       // pio::


namespace geoflow {
namespace tbox {

abort_handler EH::abort_function_ = &detail::default_abort_handler;


namespace detail {
void default_abort_handler(){
	exit(EXIT_FAILURE);
}
} // namespace detail


// STOP Nicely
void EH::abort(){
	abort_function_();
}

// Display a message and keep going
void EH::displayMessage(const std::string& Mssg, const int Line, const std::string& File){
	if( File != "" ){  pio::pout << "File: " << File << std::endl;}
	if( Line != -999 ){pio::pout << "Line: " << Line << std::endl;}
	pio::pout << Mssg << std::endl;
	pio::pout << std::endl;
}

// Display a warning message and keep going
void EH::displayWarning(const std::string& Mssg, const int Line, const std::string& File){
	pio::pout << std::endl;
	pio::pout << "*** WARNING ***" << std::endl;
	if( File != "" ){  pio::pout << "File: " << File << std::endl;}
	if( Line != -999 ){pio::pout << "Line: " << Line << std::endl;}
	pio::pout << "Message: " << Mssg << std::endl;
	pio::pout << std::endl;
}

// Display a error message and keep going
void EH::displayError(const std::string& Mssg, const int Line, const std::string& File){
	pio::perr << std::endl;
	pio::perr << "*** ERROR ***" << std::endl;
	if( File != "" ){  pio::perr << "File: " << File << std::endl;}
	if( Line != -999 ){pio::perr << "Line: " << Line << std::endl;}
	pio::perr << "Message: " << Mssg << std::endl;
	pio::perr << std::endl;
}

// Display a Bad Allocate just occurred
void EH::badNew(){
  EH::displayError("operator 'new' failed");
  EH::abort();
}

// Reset the abort call
void EH::setAbortHandler(abort_handler handle){
  abort_function_ = handle;
}

} // namespace tbox
} // namespace geoflow


