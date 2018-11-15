/*
 * range.cpp
 *
 *  Created on: Nov 14, 2018
 *      Author: bflynt
 */

#ifndef SRC_GEOFLOW_TBOX_RANGE_CPP_
#define SRC_GEOFLOW_TBOX_RANGE_CPP_


#include <iterator>    // std::iterator_traits
#include <type_traits>

//#include "geoflow/tbox/traits.hpp"

namespace geoflow {
namespace tbox {

namespace detail {

template <typename T>
struct range_iter_base {

	// Define the std::iterator traits required of all iterators
	using iterator_category = std::input_iterator_tag;
    using value_type        = T;
    using difference_type   = std::ptrdiff_t;
    using pointer           = value_type*;
    using reference         = value_type&;

    range_iter_base(value_type current = 0) : current(current) { }

    value_type operator *() const { return current; }

    value_type const* operator ->() const { return &current; }

    range_iter_base& operator ++() {
        ++current;
        return *this;
    }

    range_iter_base operator ++(int) {
        auto copy = *this;
        ++*this;
        return copy;
    }

    bool operator ==(range_iter_base const& other) const {
        return current == other.current;
    }

    bool operator !=(range_iter_base const& other) const {
        return not (*this == other);
    }

protected:
    value_type current;
};

} // namespace detail

/**
 * @brief
 * Proxy returned by range function
 *
 * This is the class returned by the range() function
 * within an range based for loop.
 */
template <typename T>
struct range_proxy {

    struct iterator : detail::range_iter_base<T> {
        iterator(T current) : detail::range_iter_base<T>(current) { }
    }; // struct iterator

    struct step_range_proxy {

        struct iter : detail::range_iter_base<T> {
            iter(T current, T step)
                : detail::range_iter_base<T>(current), step(step) { }

            using detail::range_iter_base<T>::current;

            iter& operator ++() {
                current += step;
                return *this;
            }

            iter operator ++(int) {
                auto copy = *this;
                ++*this;
                return copy;
            }

            // Loses commutativity. Iterator-based ranges are simply broken. :-(
            bool operator ==(iter const& other) const {
                return step > 0 ? current >= other.current
                                : current < other.current;
            }

            bool operator !=(iter const& other) const {
                return not (*this == other);
            }

        private:
            T step;
        }; // struct iter

        step_range_proxy(T begin, T end, T step)
            : begin_(begin, step), end_(end, step) { }

        iter begin() const { return begin_; }

        iter end() const { return end_; }

    private:
        iter begin_;
        iter end_;
    }; // struct step_range_proxy

    range_proxy(T begin, T end) : begin_(begin), end_(end) { }

    step_range_proxy step(T step) {
        return {*begin_, *end_, step};
    }

    iterator begin() const { return begin_; }

    iterator end() const { return end_; }

private:
    iterator begin_;
    iterator end_;
}; // struct range_proxy

/**
 * @brief
 * Proxy returned by infinite range function
 *
 * This is the class returned by the range() function
 * within an infinite range based for loop.
 */
template <typename T>
struct infinite_range_proxy {
    struct iterator : detail::range_iter_base<T> {
        iterator(T current = T()) : detail::range_iter_base<T>(current) { }

        bool operator ==(iterator const&) const { return false; }

        bool operator !=(iterator const&) const { return true; }
    };

    struct step_range_proxy {
        struct iter : detail::range_iter_base<T> {
            iter(T current = T(), T step = T())
                : detail::range_iter_base<T>(current), step(step) { }

            using detail::range_iter_base<T>::current;

            iter& operator ++() {
                current += step;
                return *this;
            }

            iter operator ++(int) {
                auto copy = *this;
                ++*this;
                return copy;
            }

            bool operator ==(iter const&) const { return false; }

            bool operator !=(iter const&) const { return true; }

        private:
            T step;
        };

        step_range_proxy(T begin, T step) : begin_(begin, step) { }

        iter begin() const { return begin_; }

        iter end() const { return  iter(); }

    private:
        iter begin_;
    };

    infinite_range_proxy(T begin) : begin_(begin) { }

    step_range_proxy step(T step) {
        return step_range_proxy(*begin_, step);
    }

    iterator begin() const { return begin_; }

    iterator end() const { return iterator(); }

private:
    iterator begin_;
};


/**
 * @brief
 * Range function to return indices within a range
 *
 * @details
 * Allows the usage of the range based for loop over the
 * indices of the provided range.
 *
 * \code{.cpp}
 * for (auto i : range(0,10)){
 *    cout << i << "\n";
 * }
 * \endcode
 *
 * \code{.cpp}
 * for (auto i : range(0,10).step(inc)) {
 *    cout << i << "\n";
 * }
 * \endcode
 */
template <typename T>
range_proxy<T> range(T begin, T end) {
    return {begin, end};
}

/**
 * @brief
 * Infinite range function to return indices starting at value.
 *
 * @details
 * Allows the usage of the infinite range based for loop
 * over the all indices starting at begin.
 *
 * \code{.cpp}
 * for (auto i : range(0)){
 *    cout << i << "\n";
 *    if( i == 40 ) break;
 * }
 * \endcode
 *
 * \code{.cpp}
 * for (auto i : range(0).step(inc) ){
 *    cout << i << "\n";
 *    if( i == 40 ) break;
 * }
 * \endcode
 */
template <typename T>
infinite_range_proxy<T> range(T begin) {
    return {begin};
}


} // namespace tbox
} // namespace geoflow


#endif /* SRC_GEOFLOW_TBOX_RANGE_CPP_ */
