/*
 * filter_factory.ipp
 *
 *  Created on: Nov 29, 2018
 *      Author: bryan.flynt
 */


namespace geoflow {
namespace pdeint {


template<typename ET>
typename FilterFactory<ET>::FilterBasePtr
FilterFactory<ET>::build(const tbox::PropertyTree& ptree, const std::string filtername, FilterBasePtr& filter, Grid& grid){

        // Note: the entry ptree is the main property tree!


	// Get the prop tree specified by filtername:
        PropertyTree ftree = ptree.getPropertyTree(filtername);

	// Set the default filter type
	const std::string default_filter = "none";

	// Get the type of filter
	const std::string filter_name = filtertree.getValue("filter_name", default_filter);

        // Set traits from prop tree:
        FilterTraitsType filtertraits;

        mpixx::communicator world;        
	FilterBasePtr base_ptr;
	if( "none" == filter_name ){
		using FilterImpl = NullFilter<ET>;

		// Allocate filter Implementation
		std::shared_ptr<FilterImpl> filter_impl(new FilterImpl(equation, grid));

		// Set back to base type
		base_ptr = filter_impl;
	}
        else if( "boyd_filter" == filter_name ) {
		using FilterImpl = GBoydFilter<ET>;
                typename GBoydFilter<ET>::Traits traits;

		// Allocate filter Implementation
		std::shared_ptr<FilterImpl> filter_impl(new FilterImpl(traits, grid));
		// Set back to base type
		base_ptr = filter_impl;
        }
	else {
		EH_ERROR("Requested filter not found: " << filter_name);
	}

	return base_ptr;
}



} // namespace pdeint
} // namespace geoflow

