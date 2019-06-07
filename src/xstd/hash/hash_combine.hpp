/*
 * hash_combine.hpp
 *
 *  Created on: Apr 9, 2019
 *      Author: bflynt
 */

#ifndef HASH_COMBINE_HPP_
#define HASH_COMBINE_HPP_

#include <functional>

namespace xstd{


/// Default template for combining multiple entries into a single hash
/**
 * Function that combines multiple entries into a single std::hash seed.
 * This empty function is the terminal function for the template expanded
 * function which performs the actual work.
 */
inline void hash_combine(std::size_t& seed) {
}

/// Combine multiple entries into a single hash seed value
/**
 * Function that combines multiple entries into a single std::hash seed.
 *
 * \code
 * std::size_t h=0;
 * hash_combine(h, obj1, obj2, obj3);
 * \endcode
 *
 * \note
 * The function contains magic number and bit shifts which were
 * copied from the same function within the Boost Libraries.
 */
template <typename T, typename... Rest>
inline void hash_combine(std::size_t& seed, const T& v, Rest... rest) {
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    hash_combine(seed, rest...);
}

} /* namespace xstd */

/// Overload the std::hash function to include the provided type.
/**
 * This macro overloads the std::hash<> function for a provided type
 * thereby allowing custom types to be hashed by containers which
 * require it.
 *
 * \code
 * struct MyStruct {
 *    std::string key1;
 * 	  std::string key2;
 * 	  bool key3;
 * };
 *
 * MAKE_HASHABLE(MyStruct, t.key1, t.key2, t.key3)
 *
 * std::unordered_map<MyStruct> hashed_set;
 * \endcode
 */
#ifndef MAKE_HASHABLE
#define MAKE_HASHABLE(type, ...) \
    namespace std {\
        template<> struct hash<type> {\
            std::size_t operator()(const type &t) const {\
                std::size_t ret = 0;\
                xstd::hash_combine(ret, __VA_ARGS__);\
                return ret;\
            }\
        };\
    }
#endif /* MAKE_HASHABLE */

#endif /* HASH_COMBINE_HPP_ */
