/*
 * intrusive_ptr.hpp
 *
 *  Created on: Apr 25, 2019
 *      Author: bryan.flynt
 */

#ifndef INTRUSIVE_PTR_HPP_
#define INTRUSIVE_PTR_HPP_

#include <cassert>     // assert
#include <cstddef>     // std::nullptr_t
#include <functional>  // std::less, std::hash
#include <type_traits> // std::remove_extent<T>



namespace xstd {




template<typename T>
class intrusive_ptr {

private:
	using this_type = intrusive_ptr;

public:
	using element_type = typename std::remove_extent<T>::type;


	constexpr intrusive_ptr() noexcept : data_ptr_(nullptr){}

    intrusive_ptr( T* p, bool add_ref = true ): data_ptr_(p) {
        if( (data_ptr_ != nullptr) && add_ref ){
        	intrusive_ptr_add_ref(data_ptr_);
        }
    }

    template<typename U>
    intrusive_ptr(intrusive_ptr<U> const& other): data_ptr_(other.get()){
    	if(data_ptr_ != nullptr){
    		intrusive_ptr_add_ref(data_ptr_);
    	}
    }

    intrusive_ptr(intrusive_ptr const& rhs): data_ptr_( rhs.data_ptr_ ){
        if(data_ptr_ != nullptr) intrusive_ptr_add_ref( data_ptr_ );
    }

    ~intrusive_ptr(){
        if(data_ptr_ != nullptr){
        	intrusive_ptr_release(data_ptr_);
        }
    }

    template<class U>
    intrusive_ptr& operator=(intrusive_ptr<U> const& rhs){
        this_type(rhs).swap(*this);
        return *this;
    }

    intrusive_ptr(intrusive_ptr&& rhs) noexcept : data_ptr_(rhs.data_ptr_){
        rhs.data_ptr_ = nullptr;
    }

    intrusive_ptr & operator=(intrusive_ptr&& rhs) noexcept{
        this_type(static_cast<intrusive_ptr&&>(rhs)).swap(*this);
        return *this;
    }

    template<class U> friend class intrusive_ptr;
    template<class U>
    intrusive_ptr(intrusive_ptr<U>&& rhs) noexcept : data_ptr_(rhs.data_ptr_){
    	rhs.data_ptr_ = nullptr;
    }

    template<class U>
    intrusive_ptr& operator=(intrusive_ptr<U>&& rhs) noexcept{
        this_type(static_cast<intrusive_ptr<U>&&>(rhs)).swap(*this);
        return *this;
    }

    intrusive_ptr& operator=(intrusive_ptr const& rhs){
        this_type(rhs).swap(*this);
        return *this;
    }

    intrusive_ptr& operator=(T* rhs){
        this_type(rhs).swap(*this);
        return *this;
    }

    void reset(){
        this_type().swap(*this);
    }

    void reset(T* rhs ){
        this_type(rhs).swap(*this);
    }

    void reset(T* rhs, bool add_ref ){
        this_type(rhs,add_ref).swap(*this);
    }

    T* get() const noexcept{
        return data_ptr_;
    }

    T* detach() noexcept{
        T* ret = data_ptr_;
        data_ptr_ = nullptr;
        return ret;
    }

    T& operator*() const noexcept{
        assert( data_ptr_ != nullptr );
        return *data_ptr_;
    }

    T* operator->() const noexcept{
        assert( data_ptr_ != nullptr );
        return data_ptr_;
    }

    void swap(intrusive_ptr& rhs) noexcept{
        T* tmp = data_ptr_;
        data_ptr_ = rhs.data_ptr_;
        rhs.data_ptr_ = tmp;
    }

    /// Implicit conversion to bool
	operator bool() const{
		return (data_ptr_ != nullptr);
	}

	/*
	std::size_t use_count() const noexcept{
		if(data_ptr_ != nullptr){
			return data_ptr_->use_count();
		}
		return 0;
	}
	*/

private:
	T* data_ptr_;
};


template<class T, class U> inline bool operator==(intrusive_ptr<T> const& a, intrusive_ptr<U> const &b) noexcept{
    return a.get() == b.get();
}

template<class T, class U> inline bool operator!=(intrusive_ptr<T> const& a, intrusive_ptr<U> const& b) noexcept{
    return a.get() != b.get();
}

template<class T, class U> inline bool operator==(intrusive_ptr<T> const& a, U* b) noexcept{
    return a.get() == b;
}

template<class T, class U> inline bool operator!=(intrusive_ptr<T> const& a, U* b) noexcept{
    return a.get() != b;
}

template<class T, class U> inline bool operator==(T* a, intrusive_ptr<U> const& b) noexcept{
    return a == b.get();
}

template<class T, class U> inline bool operator!=(T* a, intrusive_ptr<U> const& b) noexcept{
    return a != b.get();
}

template<class T> inline bool operator==(intrusive_ptr<T> const& p, std::nullptr_t) noexcept{
    return p.get() == nullptr;
}

template<class T> inline bool operator==(std::nullptr_t, intrusive_ptr<T> const& p) noexcept{
    return p.get() == nullptr;
}

template<class T> inline bool operator!=(intrusive_ptr<T> const& p, std::nullptr_t) noexcept{
    return p.get() != nullptr;
}

template<class T> inline bool operator!=(std::nullptr_t, intrusive_ptr<T> const& p) noexcept{
    return p.get() != nullptr;
}

template<class T> inline bool operator<(intrusive_ptr<T> const& a, intrusive_ptr<T> const& b) noexcept{
    return std::less<T*>()(a.get(), b.get());
}

template<class T> void swap(intrusive_ptr<T>& lhs, intrusive_ptr<T>& rhs) noexcept{
    lhs.swap(rhs);
}

template<class T> T* get_pointer(intrusive_ptr<T> const& p) noexcept{
    return p.get();
}

template<class T, class U> intrusive_ptr<T> static_pointer_cast(intrusive_ptr<U> const& p){
    return static_cast<T*>(p.get());
}

template<class T, class U> intrusive_ptr<T> const_pointer_cast(intrusive_ptr<U> const& p){
    return const_cast<T*>(p.get());
}

template<class T, class U> intrusive_ptr<T> dynamic_pointer_cast(intrusive_ptr<U> const& p){
    return dynamic_cast<T*>(p.get());
}

template<class Y> std::ostream & operator<<(std::ostream & os, intrusive_ptr<Y> const& p){
    os << p.get();
    return os;
}

} /* namespace xstd */



namespace std{
template<typename T>
struct hash<::xstd::intrusive_ptr<T>>{
	using argument_type = ::xstd::intrusive_ptr<T>;
	using result_type   = std::size_t;
	using key_type      = typename argument_type::element_type*;

	result_type operator()(argument_type const& ptr) const noexcept{
		return std::hash<key_type>(ptr->get());
	}
};
}



#endif /* INTRUSIVE_PTR_HPP_ */
