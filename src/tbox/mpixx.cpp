/*
 * mpixx.cpp
 *
 *  Created on: Nov 13, 2018
 *      Author: bflynt
 */

#include "tbox/mpixx.hpp"

#include <cstdlib> // exit()

namespace geoflow {
namespace tbox {
void mpixx_abort_handler(){
	mpixx::communicator world;
	world.abort(EXIT_FAILURE);
}
} // namespace tbox
} // namespace geoflow
