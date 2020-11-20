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

// Initialize the number of times the Tracer has been called
std::size_t Tracer::m_count = 0;


Tracer::Tracer(const std::string message){
#if defined( GEOFLOW_TRACER_USE_PIO )
	std::string full_text = std::string(indent(),' ') + message + " -->";
	pio::pout << full_text << std::endl;
	indent() = indent() + m_nest_indent;
#endif
#if defined( GEOFLOW_TRACER_USE_NVTX )
	m_count++;
	 constexpr std::uint32_t colors[] = {0xff00ff00, 0xff0000ff, 0xffffff00, 0xffff00ff, 0xff00ffff, 0xffff0000, 0xffffffff};
	 constexpr std::size_t num_colors = sizeof(colors)/sizeof(std::uint32_t);
	 std::size_t color_id = m_count++ %  num_colors;
	 nvtxEventAttributes_t eventAttrib = {0};
	 eventAttrib.version       = NVTX_VERSION;
     eventAttrib.size          = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
     eventAttrib.colorType     = NVTX_COLOR_ARGB;
	 eventAttrib.color         = colors[color_id];
	 eventAttrib.messageType   = NVTX_MESSAGE_TYPE_ASCII;
	 eventAttrib.message.ascii = message.c_str();
	 nvtxRangePushEx(&eventAttrib);
//	 nvtxRangePushA(message.c_str());
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
