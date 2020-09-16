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
FilterFactory<ET>::build(const tbox::PropertyTree& ptree, const std::string filterblock, Grid& grid){

        // Note: the entry ptree is the main property tree!


	// Get the prop tree specified by filterblock:
        PropertyTree ftree = ptree.getPropertyTree(filterblock);

	// Set the default filter type
	const std::string default_filter = "none";

	// Get the type of filter
	const std::string filter_name = ftree.getValue("filter_name", default_filter);


	FilterBasePtr base_ptr;
	if( "none" == filter_name 
         || ""     == filter_name ) {
		using FilterImpl = NullFilter<Types>;
#if 1
		// Allocate filter Implementation
                std::shared_ptr<FilterImpl> filter_impl(new FilterImpl());

		// Set back to base type
		base_ptr = filter_impl;
#endif
	}
        else if( "boyd_filter" == filter_name ) {
		using FilterImpl = GBoydFilter<Types>;
                typename GBoydFilter<ET>::Traits traits;

#if 0
		// Allocate filter Implementation
		std::shared_ptr<FilterImpl> filter_impl(new FilterImpl(traits, grid));
		// Set back to base type
		base_ptr = filter_impl;
#endif
        }
	else {
		EH_ERROR("Requested filter not found: " << filter_name);
	}

	return base_ptr;
}



} // namespace pdeint
} // namespace geoflow

