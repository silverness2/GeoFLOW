/*
 * property_tree.cpp
 *
 *  Created on: Nov 13, 2018
 *      Author: bflynt
 */


#include "tbox/property_tree.hpp"

#include <iostream>


using namespace geoflow::tbox;

int main() {

	PropertyTree pt;

	pt.load_file("data/generated.jsn");

	/*
	for(auto const& key : pt.getKeys() ){
		std::cout << key << std::endl;
	}

	std::cout << "Is ivalue = " << pt.isValue<int>("ivalue")     << std::endl;
	std::cout << "Is PTree  = " << pt.isValue<int>("object_1")   << std::endl;
	std::cout << "Is array  = " << pt.isValue<int>("array")      << std::endl;

	std::cout << "Is ivalue = " << pt.isArray<int>("ivalue")     << std::endl;
	std::cout << "Is PTree  = " << pt.isArray<int>("object_1")   << std::endl;
	std::cout << "Is array  = " << pt.isArray<int>("array")      << std::endl;

	std::cout << "Is ivalue = " << pt.isPropertyTree("ivalue")   << std::endl;
	std::cout << "Is PTree  = " << pt.isPropertyTree("object_1") << std::endl;
	std::cout << "Is array  = " << pt.isPropertyTree("array")    << std::endl;
	 */


	return 0;
}
