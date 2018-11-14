/*
 * tracer.hpp
 *
 *  Created on: Nov 1, 2018
 *      Author: bflynt
 */

#include "geoflow/tbox/tracer.hpp"

#if defined( MOSAIC_USE_TRACER )

#include "geoflow/tbox/pio.hpp"

namespace geoflow {
namespace tbox {

Tracer::Tracer(const std::string message){
	pio::pout << std::string(indent(),' ') << message << std::endl;
	indent() = indent() + m_nest_indent;
}

Tracer::Tracer(const std::string prefix, const std::string message){
	pio::pout << std::string(indent(),' ') << prefix << ": " << message << std::endl;
	indent() = indent() + m_nest_indent;
}

Tracer::~Tracer(){
	indent() = indent() - m_nest_indent;
	pio::pout << std::string(indent(),' ') << "---" << std::endl;
}

std::size_t& Tracer::indent(){
	static std::size_t m_current_indent{0};
	return m_current_indent;
}

} // namespace tbox
} // namespace geoflow

#endif // defined( MOSAIC_USE_TRACER )
