/*
 * morton.hpp
 *
 *  Created on: Nov 7, 2019
 *      Author: bflynt
 */

#ifndef SRC_TBOX_MORTON_HPP_
#define SRC_TBOX_MORTON_HPP_

#include <bitset>
#include <cassert>


namespace geoflow {
namespace tbox {

/** Use a std::bitset for the MortonIndex type
 *
 * We are using the standard library bitset as a
 * Morton index since it has all the features we require.
 */
template<std::size_t N>
using MortonIndex = std::bitset<N>;


/** Class to calculate the MortonIndex
 *
 * Rather than calculate the index using the full range of
 * floating point numbers this class holds a range to calculate
 * the Morton index over.
 */
template<typename COORD>
class MortonIndexer {
public:
	using Coordinate  = COORD;
	using value_type  = typename Coordinate::value_type;
	using size_type   = typename Coordinate::size_type;


	MortonIndexer() = default;
	MortonIndexer(const MortonIndexer& other) = default;
	MortonIndexer(MortonIndexer&& other) = default;
	~MortonIndexer() = default;
	MortonIndexer& operator=(const MortonIndexer& other) = default;
	MortonIndexer& operator=(MortonIndexer&& other) = default;

	MortonIndexer(const Coordinate& min_corner,
			      const Coordinate& max_corner) :
			    	  min_(min_corner), max_(max_corner) {
	}

	template<typename MortonIndex>
	void generate(const Coordinate& coords, MortonIndex& index){
		assert(min_.size() == max_.size());

		size_type  id;
		value_type mid;
		Coordinate low(min_);
		Coordinate hgh(max_);
		const size_type ndims(low.size());
		const size_type nbits(index.size());

		for(size_type i = 0; i < nbits; ++i){
			id = i % ndims;
			mid = (hgh[id] + low[id]) / value_type(2);
			if( coords[id] < mid ){
				hgh[id] = mid;
				index[i] = false;
			}
			else {
				low[id] = mid;
				index[i] = true;
			}
		}
	}

private:
	Coordinate min_;
	Coordinate max_;
};


} // namespace tbox
} // namespace geoflow


namespace std {

/** MortonIndex less than comparison
 */
template<std::size_t N>
bool operator<(const std::bitset<N>& a, const std::bitset<N>& b){
	for(std::size_t i = 0; i < a.size(); ++i){
		if(a[i] != b[i]) return b[i];
	}
	return false;
}

}

#endif /* SRC_TBOX_MORTON_HPP_ */
