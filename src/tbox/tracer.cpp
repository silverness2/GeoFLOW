/*
 * tracer.hpp
 *
 *  Created on: Nov 1, 2018
 *      Author: bflynt
 */

#include "tbox/tracer.hpp"

#if defined( GEOFLOW_USE_TRACER )


#if defined( GEOFLOW_TRACER_USE_NVTX )
#include "nvToolsExt.h"
#endif // defined( GEOFLOW_USE_NVTX )

#if defined( GEOFLOW_TRACER_USE_PIO )
#include "tbox/pio.hpp"
#endif

namespace geoflow {
namespace tbox {

Tracer::Tracer(const std::string message){
#if defined( GEOFLOW_TRACER_USE_PIO )
	std::string full_text = std::string(indent(),' ') + message + " -->";
	pio::pout << full_text << std::endl;
	indent() = indent() + m_nest_indent;
#endif
#if defined( GEOFLOW_TRACER_USE_NVTX )
	nvtxRangePushA(message.c_str());
#endif
}

Tracer::Tracer(const std::string prefix, const std::string message){
#if defined( GEOFLOW_TRACER_USE_PIO )
	std::string full_text = std::string(indent(),' ') + prefix + ": " + message;
	pio::pout << full_text << std::endl;
	indent() = indent() + m_nest_indent;
#endif
#if defined( GEOFLOW_TRACER_USE_NVTX )
	nvtxRangePushA(message.c_str());
#endif
}

Tracer::~Tracer(){
#if defined( GEOFLOW_TRACER_USE_NVTX )
	nvtxRangePop();
#endif
#if defined( GEOFLOW_TRACER_USE_PIO )
	indent() = indent() - m_nest_indent;
	pio::pout << std::string(indent(),' ') << "<--" << std::endl;
#endif
}

std::size_t& Tracer::indent(){
	static std::size_t m_current_indent{0};
	return m_current_indent;
}

} // namespace tbox
} // namespace geoflow

#endif // defined( GEOFLOW_USE_TRACER )
