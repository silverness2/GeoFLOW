/*
 * filter_factory.hpp
 *
 *  Created on: Sep. 14, 2020
 *      Author: d.rosenberg
 */

#ifndef SRC_PDEINT_FILTER_FACTORY_HPP_
#define SRC_PDEINT_FILTER_FACTORY_HPP_

#include <memory>
#include <string>
#include "tbox/error_handler.hpp"
#include "tbox/property_tree.hpp"
#include "pdeint/null_filter.hpp"
#include "pdeint/filter_base.hpp"
#include "gboyd_filter.hpp"


namespace geoflow {
namespace pdeint {


template<typename TypePack>
struct FilterFactory {
	using Types         = TypePack;
	using FilterBase    = FilterBase<Types>;
	using FilterBasePtr = std::shared_ptr<FilterBase>;
	using Grid          = typename TypePack::Grid;

	static FilterBasePtr build(const tbox::PropertyTree& ptree, std::string filtername,  FilterBasePtr& filter, Grid& grid);

};

} // namespace pdeint
} // namespace geoflow

#include "pdeint/filter_factory.ipp"

#endif /* SRC_PDEINT_FILTER_FACTORY_HPP_ */
