
#ifndef PIO_H_
#define PIO_H_

#include "configure.hpp"

#include <fstream>
#include <iostream>
#include <string>


namespace geoflow {
namespace tbox {

/**
 * Class PIO manages parallel stream I/O and logging.  Static member
 * function initialize() must be called before any of the parallel streams
 * pout, perr, or plog may be used.  Routine finalize() should also be called
 * before termination of the program.  Note that these functions should be 
 * called by the manager startup and shutdown routines and therefore
 * should not be explicitly called by an application code.
 *
 * By default, logging is disabled.  To enable logging, call one of the
 * routines logOnlyNodeZero() or logAllNodes().  Logging may be suspended
 * and resumed.
 */
struct pio{

	/**
	 * Parallel output stream pout writes to the standard output from node zero
	 * only.  Output from other nodes is ignored.  If logging is enabled, then
	 * output is mirrored to the log stream, as well.
	 */
	static std::ostream pout;

	/**
	 * Parallel output stream perr writes to the standard error from all nodes.
	 * Output is prepended with the processor number.  If logging is enabled,
	 * then output is mirrored to the log stream, as well.
	 */
	static std::ostream perr;

	/**
	 * Parallel output stream plog writes output to the log file.  When logging
	 * from multiple processors, the processor number is appended to the filename.
	 */
	static std::ostream plog;

	/**
	 * Provide an endl to pio for simplicity
	 */
	//using endl = std::endl;

   /**
    * Initialize the parallel I/O streams.  This routine must be called
    * before using pout, perr, or plog.  This routine is automatically
    * invoked by the library start-up routines.  Passing the rank into
    * the call decouples the class from the MPI implementation and the
    * requirement that MPI be initialized first.
    */
   static void initialize(std::size_t rank);

   /**
    * Shut down the parallel I/O streams and close log files.  This routine
    * must be called before program termination and is currently invoked from
    * the library shutdown procedure.
    */
   static void finalize();

   /**
    * Log messages for node zero only to the specified filename.  All output
    * to pout, perr, and plog on node zero will go to the log file.
    */
   static void logOnlyNodeZero(const std::string& filename);

   /**
    * Log messages from all nodes.  The diagnostic data for processor XXXXX
    * will be sent to a file with the name filename.XXXXX, where filename is
    * the function argument.
    */
   static void logAllNodes(const std::string& filename);

   /**
    * Temporarily suspend file logging.  Log file output will be discarded,
    * although the output file will not be closed.  Logging can be resumed
    * by calling member function resumeLogging().
    */
   static void suspendLogging();

   /**
    * Resume logging after logging was suspended via member function
    * suspendLogging().
    */
   static void resumeLogging();

  /**
   * Test if PIO has been initialized
   */
  static bool initialized(void);


private:
  static void shutdownFilestream();       // shutdown the log filestream
  static int s_rank;                      // processor rank in MPI group
  static std::ofstream* s_filestream;     // NULL or log filestream

};


} // namespace tbox
} // namespace geoflow



#endif
