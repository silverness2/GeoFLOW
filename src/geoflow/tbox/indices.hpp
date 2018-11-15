/*
 * indicies.hpp
 *
 *  Created on: Nov 14, 2018
 *      Author: bflynt
 */

#ifndef SRC_GEOFLOW_TBOX_INDICES_HPP_
#define SRC_GEOFLOW_TBOX_INDICES_HPP_

#include <iterator>    // std::iterator_traits
#include <type_traits>

#include "geoflow/tbox/range.hpp"
#include "geoflow/tbox/traits.hpp" // has_size<>

namespace geoflow {
namespace tbox {


/**
 * @brief
 * Indices of class which has a size() method
 *
 * @details
 * Allows the usage of the range based for loop over the
 * indices of the provided container with size method..
 *
 * \code{.cpp}
 * std::vector<double> vec(10);
 * for (auto i : indices(vec)){
 *    vec[i] = 3.1415 * i;
 * }
 * \endcode
 *
 * \code{.cpp}
 * std::vector<double> vec(10);
 * for (auto i : indices(vec).size(inc) ){
 *    vec[i] = 3.1415 * i;
 * }
 * \endcode
 */
template <typename C, typename = typename std::enable_if<traits::has_size<C>::value>>
auto indices(C const& cont) -> range_proxy<decltype(cont.size())> {
    return {0, cont.size()};
}


/**
 * @brief
 * Indices of simple array
 *
 * @details
 * Allows the usage of the range based for loop over the
 * indices of the provided compile time array.
 *
 * \code{.cpp}
 * double vec[3];
 * for (auto i : indices(vec)){
 *    vec[i] = 3.1415 * i;
 * }
 * \endcode
 *
 * \code{.cpp}
 * double vec[3];
 * for (auto i : indices(vec).size(inc) ){
 *    vec[i] = 3.1415 * i;
 * }
 * \endcode
 */
template <typename T, std::size_t N>
range_proxy<std::size_t> indices(T (&)[N]) {
    return {0, N};
}

/**
 * @brief
 * Indices of initializer_list
 *
 * @details
 * Allows the usage of the range based for loop over the
 * indices of the provided initializer_list
 *
 * \code{.cpp}
 * for (auto i : indices({2,3,4,5})){
 *    vec[i] = 3.1415 * i; // i = 0,1,2,3
 * }
 * \endcode
 *
 * \code{.cpp}
 * for (auto i : indices({2,3,4,5}).size(inc) ){
 *    vec[i] = 3.1415 * i; // i = 0,1,2,3
 * }
 * \endcode
 */
template <typename T>
range_proxy<typename std::initializer_list<T>::size_type>
indices(std::initializer_list<T>&& cont) {
    return {0, cont.size()};
}



} // namespace tbox
} // namespace geoflow


#endif /* SRC_GEOFLOW_TBOX_INDICES_HPP_ */
