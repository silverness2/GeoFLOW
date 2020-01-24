/*
 *io_factory.ipp
 *
 *  Created on: Jan 20, 2020 
 *      Author: d.rosenberg
 */

#include "gio.hpp"

namespace geoflow {
namespace pdeint {


template<typename ET>
typename IOFactory<ET>::IOBasePtr
IOFactory<ET>::build(const tbox::PropertyTree& ptree, Grid& grid, State& utmp){

        int itmp;

	// Set the default IO type
	const std::string default_io = "gio";

	// Get the IO object name
	std::string ioobj_name = ptree.getValue("io_implementation", default_io);

        // Set traits from prop tree:
        typename GIO<ET>::Traits gtraits;

        PropertyTree ioobj_ptree = ptree.getPropertyTree(ioobj_name);

     
	// Create the IO object and cast to base type
	IOBasePtr base_ptr;
	if( "io_method" == ioobj_name ){
		using IOImpl = GIO<ET>;

                gtraits.ivers     = ioobj_ptree.getValue <int>  ("ivers",0);
                gtraits.io_type   = ioobj_ptree.getValue <int>  ("io_type",1);
                gtraits.prgrid    = ioobj_ptree.getValue<bool>  ("prgrid",false);
                gtraits.wtime     = ioobj_ptree.getValue <int>  ("wtime",6);
                gtraits.wtask     = ioobj_ptree.getValue <int>  ("wtask",5);
                gtraits.wfile     = ioobj_ptree.getValue <int>  ("wfile",2048);
                gtraits.dim       = ioobj_ptree.getValue <int>  ("dim",GDIM);

		// Allocate IO object implementation
		std::shared_ptr<IOImpl> eqn_impl(new IOImpl(grid, gtraits, utmp));

		// Set back to base type
		base_ptr = io_impl;
	}
	else {
		EH_ERROR("Requested equation not found: " << ioobj_name);
	}

	return base_ptr;
}


} // namespace pdeint
} // namespace geoflow

