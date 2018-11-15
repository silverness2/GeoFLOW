

#ifndef GEOFLOW_ERROR_HANDLER_H_
#define GEOFLOW_ERROR_HANDLER_H_


#include <sstream>
#include <string>

#include "tbox/pio.hpp"  // pio::


namespace geoflow {
namespace tbox {


/** Definition of an abort handler function type
 */
typedef void (*abort_handler)();


namespace detail {
/** Default abort function
 *
 * This is the default function called when an abort occurs.
 * If an application requires something else such as MPI_Abort then
 * a wrapped function should be provided and passed by using
 * EH::setAbortHandler( new_abort_func );
 */
void default_abort_handler();
} // namespace detail

/**
 * @brief
 * Structure to organize display calls.
 *
 * Single structure to organize all calls for
 * display to the screen and exit if needed.
 */
struct EH{
	/**
	 * @brief
	 * Quits Program.
	 *
	 * @detail
	 * Calls appropriate stop
	 */
	static void abort();

	/**
	 * @brief
	 * Displays Message.
	 *
	 * @detail
	 * Displays a message to the screen
	 */
	static void displayMessage(const std::string& Mssg,
	                             const int Line = -999,
	                             const std::string& File = "");

	/**
	 * @brief
	 * Shows a Warning.
	 *
	 * @detail
	 * Displays a warning to the screen
	 */
	static void displayWarning(const std::string& Mssg,
	                             const int Line = -999,
	                             const std::string& File = "");

	/**
	 * @brief
	 * Shows an Error.
	 *
	 * @detail
	 * Displays an error to the screen
	 * Does not quit the program (call ABORT or EH_ERROR macro)
	 */
	static void displayError(const std::string& Mssg,
	                           const int Line = -999,
	                           const std::string& File = "");
				   
				   
	/**
	 * @brief
	 * Function to call when a bad allocate occurs
	 *
	 * @detail
	 * Function passed to the std::set_new_handler to be called
	 * when a bad new allocation occurs.  In this case it terminates
	 * the program cleanly.
	 */
  static void badNew();

  /**
   * @brief
   * Set the abort handler function.
   *
   * @detail
   * The function passed to the method will be called when ever an
   * abort is needed.  This allows decoupling of serial and MPI
   * abort.
   */
  static void setAbortHandler(abort_handler ah);

private:

  static abort_handler abort_function_;

};

} // namespace tbox
} // namespace geoflow


#define EH_MESSAGE(X)                                       \
   do {                                                     \
	  using namespace geoflow::tbox;                      \
      std::ostringstream mssg;                              \
      mssg << X << std::ends;                               \
      EH::displayMessage(mssg.str(),__LINE__,__FILE__);     \
   } while (0)


#define EH_WARNING(X)                                       \
   do {                                                     \
	  using namespace geoflow::tbox;                      \
      std::ostringstream mssg;                              \
      mssg << X << std::ends;                               \
      EH::displayWarning(mssg.str(),__LINE__,__FILE__);     \
   } while (0)


#define EH_ERROR(X)                                         \
   do {                                                     \
	  using namespace geoflow::tbox;                      \
      std::ostringstream mssg;                              \
      mssg << X << std::ends;                               \
      EH::displayError(mssg.str(),__LINE__,__FILE__);       \
      EH::abort();                                          \
   } while (0)


#endif /* GEOFLOW_ERROR_HANDLER_H_ */
