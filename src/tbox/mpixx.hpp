/*
 * mpixx.hpp
 *
 *  Created on: Nov 13, 2018
 *      Author: bflynt
 */

#ifndef SRC_GEOFLOW_MPIXX_HPP_
#define SRC_GEOFLOW_MPIXX_HPP_

#include "boost/mpi.hpp"

namespace mpixx = boost::mpi;

namespace geoflow {
namespace tbox {
void mpixx_abort_handler();
} // namespace tbox
} // namespace geoflow

#endif /* SRC_GEOFLOW_MPIXX_HPP_ */
