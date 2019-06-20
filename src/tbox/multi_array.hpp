/*
 * multi_array.hpp
 *
 *  Created on: Jun 17, 2019
 *      Author: bflynt
 */

#ifndef TEST_SCRATCH_MULTI_ARRAY_HPP_
#define TEST_SCRATCH_MULTI_ARRAY_HPP_

#include <algorithm>
#include <cassert>
#include <functional>
#include <iterator>
#include <memory>


namespace tbox {

template<typename T, std::size_t N, typename Alloc = std::allocator<T> >
class multi_array {

public:

	using value_type       = T;
	using reference        = value_type&;
	using const_reference  = const value_type&;
	using pointer          = value_type*;
	using const_pointer    = const value_type*;
	using size_type        = std::size_t;        // UnSigned/Signed Type
	using difference_type  = std::int64_t;       // Signed Type
	using index_type       = difference_type;    // Signed Type
	using allocator_type   = Alloc;
	using strider_type     = std::function<void(const size_type*,index_type*)>;

	static constexpr size_type dimensionality = N;


	multi_array();

	multi_array(const multi_array& other);

	multi_array(multi_array&& other) = default;

	~multi_array();

	multi_array& operator=(multi_array other);

	multi_array& operator=(multi_array&& other) = default;


	multi_array(const strider_type order);



	size_type num_dimension() const;

	size_type num_values() const;

	size_type size() const;

	const size_type* shape() const;

	const index_type* stride() const;

	const index_type* base() const;

	template<typename ...Args>
	reference operator()(const Args... args);

	template<typename ...Args>
	const_reference operator()(const Args... args) const;

	// Trigger a reallocation of the data
	template<typename ...Args>
	multi_array& resize(const Args... args);

	// Changes shape if new shape fits within size()
	template<typename ...Args>
	multi_array& reshape(const Args... args);

	template<typename ...Args>
	multi_array& rebase(const Args... args);

	multi_array& reorder(const strider_type order);

	void swap(multi_array& other);

	void copy(const multi_array& other);

	void fill(const value_type& value);

	void clear();

	pointer data();

	const_pointer data() const;

	pointer origin();

	const_pointer origin() const;


	static void c_order(const size_type* shape, index_type* stride){
		stride[dimensionality-1] = 1;
		if ( dimensionality > 1 ){
			for(size_type i = dimensionality-1; i --> 0;){
				stride[i] = stride[i+1] * shape[i+1];
			}
		}
	}

	static void f_order(const size_type* shape, index_type* stride){
		stride[0] = 1;
		if ( dimensionality > 1 ){
			for(size_type i = 1; i < dimensionality; ++i){
				stride[i] = stride[i-1] * shape[i-1];
			}
		}
	}

private:
	size_type                      shape_[dimensionality];
	index_type                     stride_[dimensionality];
	index_type                     base_[dimensionality];
	std::shared_ptr<value_type[]>  shared_data_ptr_;
	allocator_type                 allocator_;
	index_type                     offset_;
	strider_type                   order_;


    template<typename... Args>
    index_type index_(const Args... args) const;

    void update_stride_();

    void allocate_();

};


template<typename T, std::size_t N, typename A>
multi_array<T,N,A>::multi_array() :
	shared_data_ptr_(nullptr),
	allocator_(allocator_type()),
	offset_(0),
	order_(c_order){
	std::fill(std::begin(shape_), std::end(shape_), 0);
	std::fill(std::begin(stride_),std::end(stride_),0);
	std::fill(std::begin(base_),  std::end(base_),  0);
}

template<typename T, std::size_t N, typename A>
multi_array<T,N,A>::multi_array(const multi_array& other) :
	shared_data_ptr_(other.shared_data_ptr_),
	allocator_(other.allocator_),
	offset_(other.offset_),
	order_(c_order){
	std::copy(std::begin(other.shape_), std::end(other.shape_), std::begin(shape_));
	std::copy(std::begin(other.stride_),std::end(other.stride_),std::begin(stride_));
	std::copy(std::begin(other.base_),  std::end(other.base_),  std::begin(base_));
}

template<typename T, std::size_t N, typename A>
multi_array<T,N,A>::~multi_array(){
	this->clear();
}

template<typename T, std::size_t N, typename A>
multi_array<T,N,A>&
multi_array<T,N,A>::operator=(multi_array other){
	this->swap(other);
	return *this;
}

//template<typename T, std::size_t N, typename A>
//template<typename ...Args>
//multi_array<T,N,A>::multi_array(const Args... args) :
//	shared_data_ptr_(nullptr),
//	offset_(0),
//	order_(c_order){
//	std::fill(std::begin(base_),std::end(base_),0);
//	this->resize(args...);
//}

template<typename T, std::size_t N, typename A>
multi_array<T,N,A>::multi_array(const strider_type order) :
	shared_data_ptr_(nullptr),
	allocator_(allocator_type()),
	offset_(0),
	order_(order){
	std::fill(std::begin(shape_), std::end(shape_), 0);
	std::fill(std::begin(stride_),std::end(stride_),0);
	std::fill(std::begin(base_),  std::end(base_),  0);
}

//template<typename T, std::size_t N, typename A>
//template<typename ...Args>
//multi_array<T,N,A>::multi_array(const strider_type order, const Args... args) :
//	shared_data_ptr_(nullptr),
//	offset_(0),
//	order_(order){
//	std::fill(std::begin(base_),std::end(base_),0);
//	this->resize(args...);
//}

template<typename T, std::size_t N, typename A>
typename multi_array<T,N,A>::size_type
multi_array<T,N,A>::num_dimension() const{
	return dimensionality;
}

template<typename T, std::size_t N, typename A>
typename multi_array<T,N,A>::size_type
multi_array<T,N,A>::num_values() const{
	return this->size();
}

template<typename T, std::size_t N, typename A>
typename multi_array<T,N,A>::size_type
multi_array<T,N,A>::size() const{
	size_type ans = 1;
	for(size_type i = 0; i < dimensionality; ++i){
		ans *= shape_[i];
	}
	return ans;
}

template<typename T, std::size_t N, typename A>
const typename multi_array<T,N,A>::size_type*
multi_array<T,N,A>::shape() const{
	return this->shape_;
}

template<typename T, std::size_t N, typename A>
const typename multi_array<T,N,A>::index_type*
multi_array<T,N,A>::stride() const{
	return this->stride_;
}

template<typename T, std::size_t N, typename A>
const typename multi_array<T,N,A>::index_type*
multi_array<T,N,A>::base() const{
	return this->base_;
}

template<typename T, std::size_t N, typename A>
template<typename ...Args>
typename multi_array<T,N,A>::reference
multi_array<T,N,A>::operator()(const Args... args){
	return this->origin()[index_(args...)];
}

template<typename T, std::size_t N, typename A>
template<typename ...Args>
typename multi_array<T,N,A>::const_reference
multi_array<T,N,A>::operator()(const Args... args) const{
	return this->origin()[index_(args...)];
}

template<typename T, std::size_t N, typename A>
template<typename ...Args>
multi_array<T,N,A>&
multi_array<T,N,A>::resize(const Args... args){
	static_assert(sizeof...(args) == dimensionality);
	const size_type new_shape[] = {static_cast<size_type>(args)...};
	std::copy(std::begin(new_shape),std::end(new_shape),std::begin(shape_));
	this->update_stride_();
	this->allocate_();
	return *this;
}

template<typename T, std::size_t N, typename A>
template<typename ...Args>
multi_array<T,N,A>&
multi_array<T,N,A>::reshape(const Args... args){
	static_assert(sizeof...(args) == dimensionality);
	const size_type new_shape[] = {static_cast<size_type>(args)...};
	const auto old_size = this->size();
	std::copy(std::begin(new_shape),std::end(new_shape),std::begin(shape_));
	this->update_stride_();
	assert(old_size == this->size());
	return *this;
}

template<typename T, std::size_t N, typename A>
template<typename ...Args>
multi_array<T,N,A>&
multi_array<T,N,A>::rebase(const Args... args){
	constexpr auto arg_size = sizeof...(args);
	static_assert(arg_size == 1 || arg_size == dimensionality);
	const index_type new_base[] = {static_cast<index_type>(args)...};
	if( arg_size == 1 ){
		std::fill(std::begin(base_),std::end(base_),new_base[0]);
	}
	else {
		std::copy(std::begin(new_base),std::end(new_base),std::begin(base_));
	}
	return *this;
}

template<typename T, std::size_t N, typename A>
multi_array<T,N,A>&
multi_array<T,N,A>::reorder(const strider_type order){
	order_ = order;
	this->update_stride_();
	return *this;
}

template<typename T, std::size_t N, typename A>
void multi_array<T,N,A>::swap(multi_array& other){
	std::swap(shape_,other.shape_);
	std::swap(stride_,other.stride_);
	std::swap(base_,other.base_);
	std::swap(shared_data_ptr_,other.shared_data_ptr_);
	std::swap(allocator_,other.allocator_);
	std::swap(offset_,other.offset_);
	std::swap(order_,other.order_);
}

template<typename T, std::size_t N, typename A>
void
multi_array<T,N,A>::copy(const multi_array& other){
	std::copy(std::begin(other.shape_),  std::end(other.shape_),  std::begin(shape_));
	std::copy(std::begin(other.stride_), std::end(other.stride_), std::begin(stride_));
	std::copy(std::begin(other.base_),   std::end(other.base_),   std::begin(base_));
	allocator_ = other.allocator_;
	offset_    = other.offset_;
	order_     = other.order_;
	this->allocate_();
	const size_type sz = this->size();
	for(size_type i = 0; i < sz; ++i){
		shared_data_ptr_[i] = other.shared_data_ptr_[i];
	}
}

template<typename T, std::size_t N, typename A>
void
multi_array<T,N,A>::fill(const value_type& value){
	const size_type sz = this->size();
	for(size_type i = 0; i < sz; ++i){
		shared_data_ptr_[i] = value;
	}
}

template<typename T, std::size_t N, typename A>
void
multi_array<T,N,A>::clear() {
	shared_data_ptr_.reset();
	std::fill(std::begin(shape_),std::end(shape_),0);
	std::fill(std::begin(stride_),std::end(stride_),0);
}


template<typename T, std::size_t N, typename A>
typename multi_array<T,N,A>::pointer
multi_array<T,N,A>::data(){
	return shared_data_ptr_.get();
}

template<typename T, std::size_t N, typename A>
typename multi_array<T,N,A>::const_pointer multi_array<T,N,A>::data() const{
	return shared_data_ptr_.get();
}

template<typename T, std::size_t N, typename A>
typename multi_array<T,N,A>::pointer multi_array<T,N,A>::origin(){
	return this->data() + offset_;
}

template<typename T, std::size_t N, typename A>
typename multi_array<T,N,A>::const_pointer multi_array<T,N,A>::origin() const{
	return this->data() + offset_;
}

template<typename T, std::size_t N, typename A>
void multi_array<T,N,A>::update_stride_(){
	this->order_(shape_,stride_);
}



template<typename T, std::size_t N, typename A>
template<typename... Args>
typename multi_array<T,N,A>::index_type
multi_array<T,N,A>::index_(const Args... args) const {
	static_assert(sizeof...(args) == dimensionality);
    const index_type ijk[] = {static_cast<index_type>(args)...};
	index_type ans(0);
	for(size_type i = 0; i < dimensionality; ++i){
		assert( ijk[i] >= base_[i] );
		assert( ijk[i] <  base_[i] + shape_[i] );
		ans += stride_[i] * (ijk[i] - base_[i]);
	}
    return ans;
}

template<typename T, std::size_t N, typename A>
void
multi_array<T,N,A>::allocate_(){
	const auto new_size = this->size();
	auto dealloc = [=](T* ptr){ allocator_.deallocate(ptr,new_size); };
	shared_data_ptr_.reset(allocator_.allocate(this->size()), dealloc);
}


} // namespace tbox

#endif /* TEST_SCRATCH_MULTI_ARRAY_HPP_ */
