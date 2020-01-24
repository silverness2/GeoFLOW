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


namespace geoflow {
namespace pdeint {


template<typename IOType>
struct IOFactory {
	using IO          = IOType;
	using IOBaseType  = IOBase<Equation>;
	using IOBasePtr   = std::shared_ptr<IOBase>;
        using State       = typename IOType::State;
	using Grid        = typename IOType::Grid;

	static IOBasePtr build(const tbox::PropertyTree& ptree, Grid& grid, State& utmp);

};

} // namespace pdeint
} // namespace geoflow

#include "pdeint/io_factory.ipp"

#endif /* SRC_PDEINT_IO_FACTORY_HPP_ */
