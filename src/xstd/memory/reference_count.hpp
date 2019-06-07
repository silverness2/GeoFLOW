/*
 * intrusive.hpp
 *
 *  Created on: Apr 26, 2019
 *      Author: bflynt
 */

#ifndef REFERENCE_COUNT_HPP_
#define REFERENCE_COUNT_HPP_


#include <atomic>
#include <cassert>
#include <cstdlib>
#include <utility>


namespace xstd {

//
// Declaration of used classes and free functions
//
template <typename Base, typename Counter = std::atomic_uint32_t>
class reference_count;

template< typename Base, typename Counter>
inline void intrusive_ptr_add_ref(const reference_count<Base,Counter>* p) noexcept;

template< typename Base, typename Counter>
inline void intrusive_ptr_release(const reference_count<Base,Counter>* p) noexcept;



/// Derived class to add reference counting for usage by xstd::intrusive_ptr
/**
 * This is the derived class that all classes who intend to use
 * intrusive_ptr can create.  It adds reference counting
 * to any class and does not require reference counting to be added
 * to a class when intrusive_ptr will not be used. The inheritance
 * order makes it difficult to use for runtime polymorphic types
 * as the derived types still must inherit from a reference_count'ed
 * base.
 *
 * Example:
 * \code
 * class base {
 *    ...
 * };
 *
 * class A : public reference_count<base> {
 *    ...
 * };
 * class B : public reference_count<base> {
 *    ...
 * };
 *
 * using intrusive_base = reference_count<base>;
 *
 * xstd::intrusive_ptr<intrusive_base> vec[2];
 * vec[0] = new A();
 * vec[1] = new B();
 *
 * \endcode
 *
 * \tparam Base Class to add reference counting to
 * \tparam Counter Type of integer to use as a counter (std::atomic_uint32_t)
 */
template <typename Base, typename Counter>
class reference_count : public Base {
public:

	template<typename... Args>
	reference_count(Args&&... args) : Base(std::forward<Args>(args)...), count_(0){
	}

	~reference_count() noexcept {}

	reference_count& operator=(const reference_count&) noexcept { return *this; }

	template<typename... Args>
	reference_count& operator=(Args&&... args) noexcept {
		Base::operator=(std::forward<Args>(args)...);
		return *this;
	}

    void swap(reference_count&) noexcept {}

    std::size_t use_count() const noexcept {
        return this->count_;
    }

protected:

    friend void intrusive_ptr_add_ref<Base,Counter>(const reference_count<Base, Counter>* p) noexcept;
    friend void intrusive_ptr_release<Base,Counter>(const reference_count<Base, Counter>* p) noexcept;

private:
    mutable Counter count_;
};

/// Increment the counter by one reference
template< typename Base, typename Counter>
inline void intrusive_ptr_add_ref(const reference_count<Base,Counter>* p) noexcept{
	assert( p != nullptr );
	++(p->count_);
}

/// Decrement the counter by one reference and delete if now zero
template< typename Base, typename Counter>
inline void intrusive_ptr_release(const reference_count<Base,Counter>* p) noexcept{
	assert( p != nullptr );
	--(p->count_);
    if( p->count_ == 0 ){
        delete static_cast<Base const*>(p);
        p = nullptr;
    }
}


} /* namespace xstd */


#endif /* REFERENCE_COUNT_HPP_ */
