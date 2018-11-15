/*
 * zipp.hpp
 *
 *  Created on: Nov 14, 2018
 *      Author: bflynt
 */

#ifndef SRC_GEOFLOW_TBOX_ZIP_HPP_
#define SRC_GEOFLOW_TBOX_ZIP_HPP_


#include <iterator>    // std::iterator_traits
#include <tuple>
#include <type_traits>

#include "tbox/traits.hpp" // get_iterator<>

namespace geoflow {
namespace tbox {

namespace detail {


// Compare two tuples for item N then forward to compare N-1
template<size_t N, typename _Tuple>
struct unroll_tuple_equal {
	static bool compare( const _Tuple& a, const _Tuple& b )   {
		if ( std::get<N>( a ) == std::get<N>( b ) ) return true;
		return unroll_tuple_equal<N - 1, _Tuple>::compare( a, b );
	}
};

// Terminating comparison when N = 0
template<typename _Tuple>
struct unroll_tuple_equal<0, _Tuple>{
	static bool compare( const _Tuple& a, const _Tuple& b )   {
		return std::get<0>( a ) == std::get<0>( b );
	}
};

/**
 * @brief
 * Returns true if any elements of the tuple are equal.
 *
 * @details
 * Recursively unrolls the tuple elements to compare two at a time.
 *
 * \param first First tuple
 * \param second Second tuple
 * \returns True is there is at least one pair of elements that equal
 */
template<typename _First, typename... _Rest>
bool AnyEqualInTuple( const std::tuple<_First, _Rest...>& first, const std::tuple<_First, _Rest...>& second ){
	return unroll_tuple_equal<sizeof...( _Rest ), std::tuple<_First, _Rest...>>::compare( first, second );
}



// Call function for item N then forward to compare N-1
template<typename _Action, size_t N, typename... _TupleArgs>
struct unroll_tuple_for_each{
	static void call( std::tuple<_TupleArgs...>& tuple, _Action&& action )   {
		action.template operator()<std::tuple<_TupleArgs...>, N>( tuple );
		unroll_tuple_for_each<_Action, N - 1, _TupleArgs...>::Call( tuple, std::forward<_Action>( action ) );
	}
};

// Terminating function when N = 0
template<typename _Action, typename... _TupleArgs>
struct unroll_tuple_for_each<_Action, 0, _TupleArgs...>{
	static void call( std::tuple<_TupleArgs...>& tuple, _Action&& action )   {
		action.template operator()<std::tuple<_TupleArgs...>, 0>( tuple );
	}
};


/**
 * @brief
 * Calls a function for each element in the given tuple.
 *
 * \param tuple The tuple on which to call the function
 * \tparam _Action Function template with operator() to take a parameter of each of the tuple types
 * \tparam _FuncArgs Additional arguments for the _Action functor
 * \tparam _FirstArg First argument of the tuple
 * \tparam _Rest Optional other tuple arguments
 */
template<typename _Action, typename... _TupleArgs>
void ForEachInTuple( std::tuple<_TupleArgs...>& tuple, _Action&& action ){
	//TODO Explicitly forbid empty tuples
	static_assert( sizeof...( _TupleArgs ) > 0, "Empty tuple is not allowed!" );
	unroll_tuple_for_each<_Action, sizeof...( _TupleArgs ) - 1, _TupleArgs...>::call( tuple, std::forward<_Action>( action ) );
}

/**
 * @brief
 * Generate a sequence of indices at compile time.
 */
template<int... S>
struct Sequence{};

template<int N, int... S>
struct SequenceGenerator : SequenceGenerator<N - 1, N - 1, S...> {};

template<int... S>
struct SequenceGenerator<0, S...>{
	using type = Sequence<S...>;
};

/**
 * @brief
 * Pass-through function used to expand a function call on all arguments of a parameter pack
 */
template<typename... T>
inline void PassThrough( T&&... ) {}

} // namespace detail

/**
 * @brief
 * Class to hold tuple of iterators with ability to increment them all.
 *
 */
template<typename... IterTypes>
struct iterator_collection {

	template<size_t Index>
	using value_type_t = typename std::tuple_element<Index, std::tuple<IterTypes...>>::type::value_type;

	using value_ref_tuple_t = std::tuple<typename std::iterator_traits<IterTypes>::reference...>;

	iterator_collection( IterTypes&&... iterators ) :
		iterator_pack_( std::forward<IterTypes>(iterators)... ){
	}

	inline bool MatchAny( const iterator_collection& other ) const      {
		return detail::AnyEqualInTuple( iterator_pack_, other.iterator_pack_ );
	}

	inline bool operator==( const iterator_collection& other ) const
    		  {
		return iterator_pack_ == other.iterator_pack_;
    		  }

	inline bool operator!=( const iterator_collection& other ) const
    		  {
		return !operator==( other );
    		  }

	inline value_ref_tuple_t Deref()
	{
		return DerefInternal( typename detail::SequenceGenerator<sizeof...( IterTypes )>::type() );
	}

	inline void Increment()
	{
		IncrementInternal( typename detail::SequenceGenerator<sizeof...( IterTypes )>::type() );
	}

private:
	template<int... S>
	inline value_ref_tuple_t DerefInternal( detail::Sequence<S...> )
	{
		return value_ref_tuple_t( *std::get<S>( iterator_pack_ )... );
	}

	template<int... S>
	inline void IncrementInternal( detail::Sequence<S...> )
	{
		detail::PassThrough( std::get<S>( iterator_pack_ ).operator++( )... );
	}

	std::tuple<IterTypes...> iterator_pack_;
};


/**
 * @brief
 * Zip Iterator which iterates over all iterators in the tuple
 *
 * @details
 * Returns a tuple of elements at the current position when dereferenced. Since
 * the collections might be of different lengths, this iterator stops when the
 * collection with the fewest elements is exhausted.
 */
template<typename... _Iters>
class zip_iterator {
	using iter_collection_type = iterator_collection<_Iters...>;
	using value_ref_tuple_t = std::tuple<typename std::iterator_traits<_Iters>::reference...>;
public:
	zip_iterator( iter_collection_type cur ) :
		current_( cur ){
}

	inline value_ref_tuple_t operator*( )
	{
		return current_.Deref();
	}

	inline zip_iterator& operator++( )
    		  {
		current_.Increment();
		return *this;
    		  }

	inline bool operator==( const zip_iterator& other ) const
    		  {
		//Again, for the comparison inside a range based for loop, one match is enough!
		return current_.MatchAny( other.current_ );
    		  }

	inline bool operator!=( const zip_iterator& other ) const
    		  {
		return !operator==( other );
    		  }

private:
	iter_collection_type current_;
};

/**
 * @brief
 * Proxy class that zips multiple iterators together
 *
 * @details
 *
 */
template<typename... _Iters>
class zip_proxy {

	using iter_collection_type = iterator_collection<_Iters...>;

public:
	zip_proxy( iter_collection_type&& begins, iter_collection_type&& ends ) :
		begin_( std::forward<iter_collection_type>( begins ) ),
		end_( std::forward<iter_collection_type>( ends ) ){
	}

	inline zip_iterator<_Iters...> begin( )   {
		return zip_iterator<_Iters...>( begin_ );
	}

	inline zip_iterator<_Iters...> end( )      {
		return zip_iterator<_Iters...>(end_ );
	}

private:
	iter_collection_type begin_;
	iter_collection_type end_;
};

/**
 * @brief
 * Creates a zip iterator to iterator over a range of collections simultaneously
 *
 * \param args All the collections to iterator over
 * \tparam Args Types of collections
 * \returns A zip_iterator over all the collections
 *
 * \code{.cpp}
 * std::vector<atype> a(sz);
 * std::list<btype>   b(sz);
 * for(auto i : zip(a,b)){
 * 	atype val_a = std::get<0>(i);
 * 	btype val_b = std::get<1>(i); *
 * }
 * \endcode
 *
 * \code{.cpp}
 * std::vector<atype>     a(100);
 * std::list<btype>       b(100);
 * std::array<ctype,100>  c(100);
 * for(auto i : zip(a,b,c).step(inc)){
 * 	atype val_a = std::get<0>(i);
 * 	btype val_b = std::get<1>(i);
 * 	ctype val_c = std::get<2>(i);
 * }
 * \endcode
 *
 */
template<typename... Args>
zip_proxy<typename traits::get_iterator<Args>::type...> zip( Args&&... args ) {
	using iter_collection_type = iterator_collection<typename traits::get_iterator<Args>::type...>;
	return zip_proxy<typename traits::get_iterator<Args>::type...>( iter_collection_type( std::begin( args )... ),
			iter_collection_type( std::end(args)... ) );
}


} // namespace tbox
} // namespace geoflow



#endif /* SRC_GEOFLOW_TBOX_ZIP_HPP_ */
