/**
 * \file traits.hpp
 * \brief Defines helpful trait tests for common tasks.
 *
 * Common tests which are not included in the std::type_traits header.want to
 *
 * \author Bryan Flynt
 */
#ifndef STDX_TRAITS_HPP_
#define STDX_TRAITS_HPP_


#include <type_traits>


namespace xstd {

/// Test if class has size() method
/**
 * Test if template type C has size() method.
 *
 * Usage:
 * \code
 * has_size<MyClass>::type
 * has_size<MyClass>::value
 * \endcode
 *
 * \tparam C Class to check
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


/// Test if class has iterator type
/**
 * Test if template type C has an associated iterator type.
 *
 * Usage:
 * \code
 * has_iterator<MyClass>::type
 * has_iterator<MyClass>::value
 * \endcode
 *
 * \tparam C Class to check
 */
template <typename C>
struct has_iterator {
    template <typename T>
    static auto check(T*) -> std::is_same<decltype(T::iterator),decltype(T::iterator)>;

    template <typename>
    static auto check(...) -> std::false_type;

    using type = decltype(check<C>(0));
    static constexpr bool value = type::value;
};

/// Test if class has begin() method
/**
 * Test if template type C has a begin() method.
 *
 * Usage:
 * \code
 * has_begin<MyClass>::type
 * has_begin<MyClass>::value
 * \endcode
 *
 * \tparam C Class to check
 */
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

/// Test if class has end() method
/**
 * Test if template type C has a end() method.
 *
 * Usage:
 * \code
 * has_end<MyClass>::type
 * has_end<MyClass>::value
 * \endcode
 *
 * \tparam C Class to check
 */
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

/// Test if class has begin() and end() methods
/**
 * Test if template type C has both begin() and end() methods.
 *
 * Usage:
 * \code
 * has_begin_end<MyClass>::type
 * has_begin_end<MyClass>::value
 * \endcode
 *
 * \tparam C Class to check
 */
template <typename C>
struct has_begin_end {
    using type = typename std::conditional<has_begin<C>::value,
    		typename has_end<C>::type, std::false_type>::type;
    static constexpr bool value = type::value;
};

/// Get the proper constantness iterator for the class
/**
 * Return the proper iterator for a class instance based
 * on the constantness of the template type.  Used when a function,
 * etc. does not know if the variable passed in a constant or not.
 *
 * Results:
 * | Input Type | Result Type    |
 * | :--------: | :------------: |
 * | T          | iterator       |
 * | const T    | const_iterator |
 *
 * Usage:
 * \code
 *	template<typename T>
 *	int do_something(T value){
 *		get_iterator<T>::type it   = T.begin();
 *		get_iterator<T>::type last = T.end();
 *		while( it != last ){
 *		   ...
 *		   ++it;
 *		}
 *	};
 * \endcode
 *
 * \tparam T Class to check
 */
template<typename T>
struct get_iterator {
  using type = typename std::conditional< std::is_const<T>::value,
      typename std::remove_reference<T>::type::const_iterator,
      typename std::remove_reference<T>::type::iterator >::type;
};

} /* namespace xstd */



#endif /* STDX_TRAITS_HPP_ */
