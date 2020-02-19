/*
 * observer_factory.hpp
 *
 *  Created on: Nov 29, 2018
 *      Author: bryan.flynt
 */

#ifndef SRC_PDEINT_OBSERVER_FACTORY_HPP_
#define SRC_PDEINT_OBSERVER_FACTORY_HPP_

#include <memory>
#include <string>
#include "tbox/error_handler.hpp"
#include "pdeint/observer_base.hpp"
#include "tbox/property_tree.hpp"
#include "pdeint/null_observer.hpp"
#include "pdeint/io_base.hpp"
#include "gio_observer.hpp"
#include "io_factory.hpp"
#include "gglobaldiag_basic.hpp"


namespace geoflow {
namespace pdeint {


template<typename EquationType>
struct ObserverFactory {
	using Equation    = EquationType;
        using EqnBase     = EquationBase<EquationType>;
        using EqnBasePtr  = std::shared_ptr<EqnBase>;
        using IOBaseType  = IOBase<EquationType>;
        using IOBasePtr   = std::shared_ptr<IOBaseType>;
	using ObsBase     = ObserverBase<Equation>;
	using ObsBasePtr  = std::shared_ptr<ObsBase>;
	using Grid        = typename EquationType::Grid;

	static ObsBasePtr build(const tbox::PropertyTree& ptree, std::string obsname,  EqnBasePtr& equation, Grid& grid, IOBasePtr pIO);

};

} // namespace pdeint
} // namespace geoflow

#include "pdeint/observer_factory.ipp"

#endif /* SRC_PDEINT_OBSERVER_FACTORY_HPP_ */
