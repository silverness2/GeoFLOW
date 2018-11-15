/*
 * enumerate.hpp
 *
 *  Created on: Nov 14, 2018
 *      Author: bflynt
 */

#ifndef SRC_GEOFLOW_TBOX_ENUMERATE_HPP_
#define SRC_GEOFLOW_TBOX_ENUMERATE_HPP_

#include <iterator> // std::iterator_traits
#include <utility>  // std::pair<>

#include "tbox/traits.hpp" // has_begin_end<>


namespace geoflow {
namespace tbox {


/**
 * @brief
 * Proxy class that enumerates iterators
 *
 * @details
 *
 */
template<typename IterType>
struct enumerate_proxy {
	using difference_type = typename std::iterator_traits<IterType>::difference_type;

	struct iterator {
		using iter_type            = IterType;
		using iter_value_type      = typename std::iterator_traits<iter_type>::value_type;
		using iter_reference_type  = typename std::iterator_traits<iter_type>::reference;

		using iterator_category   = typename std::iterator_traits<iter_type>::iterator_category;
		using difference_type     = typename std::iterator_traits<iter_type>::difference_type;
		using value_type          = std::pair<const difference_type, iter_reference_type>;
		using pointer             = value_type*;
		using reference           = value_type&;


		iterator(IterType current, difference_type index) :
			current_(current), index_(index) {
		}

		iterator& operator++() {
			++current_;
			++index_;
			return *this;
		}

		iterator operator ++(int) {
			auto copy = *this;
			++(*this);
			return copy;
		}

		bool operator ==(iterator const& other) const {
			return current_ == other.current_;
		}

		bool operator !=(iterator const& other) const {
			return not (*this == other);
		}

		bool operator <(iterator const& other) const {
			return (current_ < other.current_);
		}


		value_type operator*() { return {index_, *current_}; }

		//value_type const* operator ->() const { return &current_; }

	private:
		iter_type       current_;
		difference_type index_;
	};

	struct step_enumerate_proxy {
		struct step_iterator {
			using iter_type           = iterator;
			using iter_value_type     = typename std::iterator_traits<iter_type>::value_type;
			using iter_reference_type = typename std::iterator_traits<iter_type>::reference;

			using iterator_category   = typename std::iterator_traits<iter_type>::iterator_category;
			using difference_type     = typename std::iterator_traits<iter_type>::difference_type;
			using value_type          = iter_value_type;
			using pointer             = value_type*;
			using reference           = value_type&;


			step_iterator(iterator current, difference_type step) :
				current_(current), step_(step) {
			}

			step_iterator& operator++() {
				for(difference_type i = 0; i < step_; ++i){
					++current_;
				}
				return *this;
			}

			step_iterator operator ++(int) {
				auto copy = *this;
				++(*this);
				return copy;
			}

			bool operator ==(step_iterator const& other) const {
				 // Protect against skipping over other
				return step_ > 0 ? (other.current_ < current_ || other.current_ == current_)
						         : current_ < other.current_;
			}

			bool operator !=(step_iterator const& other) const {
				return not (*this == other);
			}

			value_type operator*() { return *current_; }

			//value_type const* operator ->() const { return &current_; }

		private:
			iter_type             current_;
			const difference_type step_;
		};

		step_enumerate_proxy(iterator begin, iterator end, difference_type step) :
			begin_(begin,step),
			end_(end,step) {
		}

		step_iterator begin() const { return begin_; }

		step_iterator end() const { return end_; }

	private:
		step_iterator        begin_;
		step_iterator        end_;
	};


	enumerate_proxy(IterType begin, IterType end, difference_type initial) :
		begin_(begin,initial),
		end_(end,end-begin) {
	}

	step_enumerate_proxy step(difference_type step) {
		return {begin_, end_, step};
	}

	iterator begin() const { return begin_; }

	iterator end() const { return end_; }

private:
	iterator        begin_;
	iterator        end_;
};


/**
 * @brief
 * Creates a enumerated iterator to iterator from first to last
 *
 *
 * \param first Starting iterator
 * \param last  End iterator which will not be included
 * \param initial Starting offset to start iterating from
 *
 * \tparam Iterator The iterator type of the container
 *
 * \returns A enumerate_proxy over all the iterators
 *
 * \code{.cpp}
 * std::vector<dtype> pvec(sz);
 * for (auto i : enumerate(pvec.begin(),pvec.end()) ){
 * 		int   index = std::get<0>(i);
 * 		dtype myval = std::get<1>(i);
 * }
 * \endcode
 *
 * \code{.cpp}
 * std::vector<dtype> pvec(sz);
 * for (auto i : enumerate(pvec.begin(),pvec.end()).step(inc) ){
 * 		int   index = std::get<0>(i);
 * 		dtype myval = std::get<1>(i);
 * }
 * \endcode
 */
template<typename Iterator>
enumerate_proxy<Iterator> enumerate( Iterator first,
		                             Iterator last,
		                             std::ptrdiff_t initial = 0) {
	return {first, last, initial};
}

/**
 * @brief
 * Creates an enumerated iterator over all elements of type
 * with begin and end methods.
 *
 * \param content Container instance to iterate over
 *
 * \tparam Container type to iterate over
 *
 * \returns A enumerate_proxy over all the iterators
 *
 * \code{.cpp}
 * std::vector<dtype> pvec(sz);
 * for (auto i : enumerate(pvec) ){
 * 		int   index = std::get<0>(i);
 * 		dtype myval = std::get<1>(i);
 * }
 * \endcode
 *
 * \code{.cpp}
 * std::vector<dtype> pvec(sz);
 * for (auto i : enumerate(pvec).step(inc) ){
 * 		int   index = std::get<0>(i);
 * 		dtype myval = std::get<1>(i);
 * }
 * \endcode
 */
template<typename Container, typename = typename std::enable_if<traits::has_begin_end<Container>::value> >
auto enumerate(Container &content) -> enumerate_proxy< decltype(content.begin()) > {
	return {content.begin(), content.end(), 0};
}

/**
 * @brief
 * Creates an enumerated iterator over all elements of an initializer_list.
 *
 * \param content Container type to iterate over
 *
 * \tparam T Type contained within initializer_list
 *
 * \returns A enumerate_proxy over all the iterators
 *
 * \code{.cpp}
 * for (auto i : enumerate( {val1, val2, val3} ) ){
 * 		int   index = std::get<0>(i);
 * 		dtype myval = std::get<1>(i);
 * }
 * \endcode
 *
 * \code{.cpp}
 * for (auto i : enumerate( {val1, val2, val3} ).step(inc) ){
 * 		int   index = std::get<0>(i);
 * 		dtype myval = std::get<1>(i);
 * }
 * \endcode
 */
template<typename T>
auto enumerate(std::initializer_list<T>&& content) -> enumerate_proxy< decltype(std::begin(content)) >{
	return {std::begin(content), std::end(content), 0};
}

/**
 * @brief
 * Creates an enumerated iterator over all elements of a static array.
 *
 * \param x Static array to iterate over
 *
 * \tparam T Type contained within static array
 *
 * \returns A enumerate_proxy over all the iterators
 *
 * \code{.cpp}
 * dtype vec[100];
 * for (auto i : enumerate(vec) ){
 * 		int   index = std::get<0>(i);
 * 		dtype myval = std::get<1>(i);
 * }
 * \endcode
 *
 * \code{.cpp}
 * dtype vec[100];
 * for (auto i : enumerate(vec).step(inc) ){
 * 		int   index = std::get<0>(i);
 * 		dtype myval = std::get<1>(i);
 * }
 * \endcode
 */
template <typename T, std::size_t N>
enumerate_proxy<T*> enumerate(T (&x)[N]){
	return {&(x[0]), &(x[0])+N, 0};
}


} // namespace tbox
} // namespace geoflow

#endif /* SRC_GEOFLOW_TBOX_ENUMERATE_HPP_ */
