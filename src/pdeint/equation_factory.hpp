/*
 * equation_factory.hpp
 *
 *  Created on: July 8, 2019 
 *      Author: bryan.flynt, d.rosenberg
 */

#ifndef SRC_PDEINT_EQUATION_FACTORY_HPP_
#define SRC_PDEINT_EQUATION_FACTORY_HPP_

#include <memory>
#include <string>
#include "gburgers.hpp"
#include "gmconv.hpp"
#include "tbox/error_handler.hpp"
#include "pdeint/equation_base.hpp"
#include "pdeint/filter_factory.hpp"
#include "tbox/property_tree.hpp"


namespace geoflow {
namespace pdeint {


template<typename EquationType>
struct EquationFactory {
	using Equation    = EquationType;
	using EqnBase     = EquationBase<Equation>;
	using EqnBasePtr  = std::shared_ptr<EqnBase>;
        using State       = typename EquationType::State;
	using Grid        = typename EquationType::Grid;
        using FilterBasePtr = std::shared_ptr<FilterBase<EquationType>>;
        using FilterList    = std::vector<FilterBasePtr>;


	static EqnBasePtr build(const tbox::PropertyTree& ptree, Grid& grid);
        static void       config_filters(const PropertyTree& ptree, const std::string eq_name,  Grid& grid, FilterList& filter_list);

};

} // namespace pdeint
} // namespace geoflow

#include "pdeint/equation_factory.ipp"

#endif /* SRC_PDEINT_EQUATION_FACTORY_HPP_ */
