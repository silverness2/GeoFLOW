/*
 * global_manager.hpp
 *
 *  Created on: Nov 14, 2018
 *      Author: bflynt
 */

#ifndef SRC_GEOFLOW_TBOX_GLOBAL_MANAGER_HPP_
#define SRC_GEOFLOW_TBOX_GLOBAL_MANAGER_HPP_


namespace geoflow {
namespace tbox {

/**
 * @brief Utility for managing startup and shutdown of all objects.
 *
 * The startup/shutdown mechanism in GlobalManager is used to manage the
 * allocation and deallocation of certain memory, particularly static data
 * that must exist during the full extent of a run or for the full extent
 * of a single problem within a run.  The specific data that is controlled
 * by this mechanism is managed using the StartupShutdownManager.
 *
 * The four steps of the startup/shutdown mechanism are:
 *
 * <ul> initialize -- called at the start of a program after MPI is
 *                    initialized but before any other objects are used.
 * <ul> startup -- called to begin a problem-specific segment of the code.
 * <ul> shutdown -- called at the end of a problem-specific segment of the
 *                  code.  Shuts down and deallocates everything that was
 *                  started and allocated by startup.
 * <ul> finalize -- called at the end of a program right before MPI is
 *                  finalized.
 *
 * The startup and shutdown functions may be called multiple times within a
 * run, in order to allow for the execution of more than one problem within one
 * program.  Initialize and finalize must be called exactly once in a single
 * run.
 *
 * @see StartupShutdownManager
 */
class GlobalManager{
public:
	/**
	 * @brief Initial setup of the package.
	 *
	 * This function should be invoked ONLY ONCE at the start of a process
	 * to initialize and AFTER MPI is initialized (if used) by a
	 * call to one of the MPI init routines.
	 *
	 * This function initializes I/O, as well as data for any classes
	 * that implement the initialize call back
	 * interface through StartupShutdownManager.
	 */
	static void initialize(int argc, char* argv[]);

	/**
	 * @brief Startup of the package.
	 *
	 * This function invokes startup for any classes that implement the
	 * startup callback interface through StartupShutdownManager.
	 * This function may be invoked more than once in a process if
	 * solving multiple problems.
	 */
	static void startup();

	/**
	 * @brief Shutdown the package.
	 *
	 * This function invokes shutdown for any classes that implement the
	 * startup callback interface through StartupShutdownManager.
	 * This function may be invoked more than once in an process if
	 * solving multiple problems.
	 */
	static void shutdown();

	/**
	 * @brief Final cleanup of the package.
	 *
	 * This function should be invoked ONLY ONCE at the end of a process
	 * to complete the cleanup of memory allocations and
	 * any other cleanup tasks. I/O will be finalized, as well as
	 * data for any classes that implement the finalize callback
	 * interface through StartupShutdownManager.
	 *
	 * After this function is called, the only thing that should occur before
	 * exiting the program is a call to MPI finalize().
	 *
	 * This function should be invoked only once.
	 */
	static void finalize();

	/**
	 * @brief Returns true if GlobalManager has been initialized.
	 */
	static bool isInitialized();

	/**
	 * @brief Returns true if GlobalManager has been started.
	 */
	static bool isStarted();

private:

   /**
    * Flag indicating GlobalManager has been initialized.
    */
   static bool s_initialized;

   /**
    * Flag indicating startup has occurred.
    */
   static bool s_started;

};

} // namespace tbox
} // namespace geoflow

#endif /* SRC_GEOFLOW_TBOX_GLOBAL_MANAGER_HPP_ */
