/*
 * lin_solver_base.ipp
 *
 *  Created on: Mar. 16, 2020
 *      Author: D. Rosenberg
 */


namespace geoflow {
namespace pdeint {


template<typename Types>
typename LinSolverBase<Types>::GNormType
LinSolverBase<Types>::str2normtype(std::string snorm)
{
  std::string s0;
  for ( auto j=0; j<GCG_NORM_NONE-1; j++ ) {
    s0 = sGBdyType[j];
    if ( snorm.compare(s0) == 0 ) return static_cast<GNormType>(j);
  }
  assert(FALSE && "Invalid norm type");

}


} // namespace pdeint
} // namespace geoflow

