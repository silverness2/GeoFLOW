/*
 * traits.hpp
 *
 *  Created on: Nov 14, 2018
 *      Author: bflynt
 */

#ifndef SRC_GEOFLOW_TBOX_TRAITS_HPP_
#define SRC_GEOFLOW_TBOX_TRAITS_HPP_



#include <type_traits>

namespace geoflow {
namespace tbox {
namespace traits {


/**
 * @brief
 * Test if template type C has size() method.
 *
 */
template <typename C>
struct has_size {
    template <typename T>
    static auto check(T*) ->
		typename std::is_integral<decltype(std::declval<T const>().size())>::type;

    template <typename>
    static auto check(...) -> std::false_type;

    using type = decltype(check<C>(0));
    static constexpr bool value = type::value;
};

template <typename C>
struct has_iterator {
    template <typename T>
    static auto check(T*) -> std::is_same<decltype(T::iterator),decltype(T::iterator)>;

    template <typename>
    static auto check(...) -> std::false_type;

    using type = decltype(check<C>(0));
    static constexpr bool value = type::value;
};

template <typename C>
struct has_begin {
    template <typename T>
    static auto check(T*) -> std::is_same<decltype(std::declval<T>().begin()),
	                                      decltype(std::declval<T>().begin())>;

    template <typename>
    static auto check(...) -> std::false_type;

    using type = decltype(check<C>(0));
    static constexpr bool value = type::value;
};

template <typename C>
struct has_end {
    template <typename T>
    static auto check(T*) -> std::is_same<decltype(std::declval<T>().end()),
	                                      decltype(std::declval<T>().end())>;

    template <typename>
    static auto check(...) -> std::false_type;

    using type = decltype(check<C>(0));
    static constexpr bool value = type::value;
};

template <typename C>
struct has_begin_end {
    using type = typename std::conditional<has_begin<C>::value,
    		typename has_end<C>::type, std::false_type>::type;
    static constexpr bool value = type::value;
};

// Take type T and return Iterator for the type.
// If input Type       => return T::iterator
// If input const Type => return T::const_iterator
template<typename T>
struct get_iterator {
  using type = typename std::conditional< std::is_const<T>::value,
      typename std::remove_reference<T>::type::const_iterator,
      typename std::remove_reference<T>::type::iterator >::type;
};

} // namespace traits
} // namespace tbox
} // namespace geoflow



#endif /* SRC_GEOFLOW_TBOX_TRAITS_HPP_ */
