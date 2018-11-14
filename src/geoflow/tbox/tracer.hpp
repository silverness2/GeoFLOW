/*
 * tracer.hpp
 *
 *  Created on: Dec 23, 2017
 *      Author: bflynt
 */

#ifndef MOSAIC_DEBUG_TRACER_HPP_
#define MOSAIC_DEBUG_TRACER_HPP_

#include "geoflow/configure.hpp"

/** \file tracer.hpp
 * Implementation of Tracer tools to trace program execution.
 *
 * The Tracer tools can be used to trace program execution by
 * printing to screen a message indented by the nested creation
 * of Tracer classes.  To enable zero runtime cost for release
 * builds a macro is used to only insert Tracer construction
 * when TRACER_ON is defined.
 */
#if defined( MOSAIC_USE_TRACER )

#include <string>


namespace geoflow {
namespace tbox {

/**
 * \brief Tracer to track program execution to stdout
 *
 * Tracer will trace a programs execution to the screen
 * when a Tracer is constructed and allowed to be destructed
 * when it goes out of scope.
 *
 * Example:
 * \code{.cpp}
 * int MyFunction(){
 *    Tracer my_local_tracer("My Message");
 *
 *    // Program Contents
 * };
 * \endcode
 *
 */
class Tracer {

public:

	Tracer() = delete;
	Tracer(const Tracer& T) = delete;
	Tracer& operator=(const Tracer& T) = delete;

	/**
	 * \brief
	 * Construct with a message and no prefix.
	 */
	Tracer(const std::string message);

	/**
	 * \brief
	 * Construct with a prefix and message.
	 */
	Tracer(const std::string prefix, const std::string message);

	/**
	 * \brief
	 * Destructor will decrease indention
	 */
	~Tracer();

	/**
	 * \brief
	 * Returns a reference to the current line indention
	 * count.
	 *
	 * \details
	 * Returns a reference to the size of the current
	 * line indention.  This is determined by the current
	 * level of nesting.  This is implemented as a static
	 * function since static data members cannot be
	 * initialized within header files.
	 */
	static std::size_t& indent();

private:
	static constexpr std::size_t m_nest_indent = 3;
	static std::size_t m_current_indent;
};

} // namespace tbox
} // namespace geoflow



/**
 * \def TRACER(msg)
 * \brief A macro to insert a tracer with message
 * \param	msg	Message to insert after prefix.
 */

/**
 * \def TRACER_W_PREFIX(pre,msg)
 * \brief A macro to insert a tracer with prefix and message
 * \param	pre	Prefix to insert before message.
 * \param	msg	Message to insert after prefix.
 */
#ifndef TRACER
   #define TRACER(msg) geoflow::tbox::Tracer local_scope_tracer(msg)
   #define TRACER_W_PREFIX(pre,msg) geoflow::tbox::Tracer local_scope_tracer(pre,msg)
#endif

#else // #ifdef MOSAIC_USE_TRACER

/**
 * Macro to ignore Tracer construction
 */
   #ifndef TRACER
      #define TRACER(msg)
      #define TRACER_W_PREFIX(pre,msg)
   #endif

#endif //  MOSAIC_USE_TRACER


#endif /* MOSAIC_DEBUG_TRACER_HPP_ */
