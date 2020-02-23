/*
 *io_factory.hpp
 *
 *  Created on: Jan 20, 2020 
 *      Author: d.rosenberg
 */

#ifndef SRC_PDEINT_IO_FACTORY_HPP_
#define SRC_PDEINT_IO_FACTORY_HPP_

#include <memory>
#include <string>
#include "tbox/error_handler.hpp"
#include "pdeint/io_base.hpp"
#include "tbox/property_tree.hpp"
#include "gio.hpp"


namespace geoflow {
namespace pdeint {


template<typename TypePack>
struct IOFactory {
	using IO          = TypePack;
	using IOBaseType  = IOBase<IO>;
	using IOBasePtr   = std::shared_ptr<IOBaseType>;
        using State       = typename TypePack::State;
	using Grid        = typename TypePack::Grid;

	static IOBasePtr build(const tbox::PropertyTree& ptree, Grid& grid, GC_COMM comm);

};

} // namespace pdeint
} // namespace geoflow

#include "pdeint/io_factory.ipp"

#endif /* SRC_PDEINT_IO_FACTORY_HPP_ */
