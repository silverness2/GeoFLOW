/*
 * array_math.hpp
 *
 *  Created on: May 21, 2019
 *      Author: bflynt
 */

#ifndef ARRAY_MATH_HPP_
#define ARRAY_MATH_HPP_


#include <array>
#include <cstdlib>
#include <type_traits>


// ============================================================
//                    Unary Operations
// ============================================================
template<typename T1, std::size_t N1>
std::array<T1,N1> operator -(const std::array<T1,N1>& a) noexcept{
	auto ans(a);
	for(std::size_t i = 0; i < N1; ++i){
		a[i] = -a[i];
	}
	return ans;
}
template<typename T1, std::size_t N1>
std::array<T1,N1> operator +(const std::array<T1,N1>& a) noexcept{
	return a;
}

// ============================================================
//                Array / Scalar Operations
// ============================================================
template<typename T1, std::size_t N1, typename T2>
void operator +=(std::array<T1,N1>& a, const T2& b) noexcept{
	static_assert(std::is_convertible<T2,T1>::value, "Type Mismatch");
	for(std::size_t i = 0; i < N1; ++i){
		a[i] += b;
	}
}

template<typename T1, std::size_t N1, typename T2>
void operator -=(std::array<T1,N1>& a, const T2& b) noexcept{
	static_assert(std::is_convertible<T2,T1>::value, "Type Mismatch");
	for(std::size_t i = 0; i < N1; ++i){
		a[i] -= b;
	}
}

template<typename T1, std::size_t N1, typename T2>
void operator *=(std::array<T1,N1>& a, const T2& b) noexcept{
	static_assert(std::is_convertible<T2,T1>::value, "Type Mismatch");
	for(std::size_t i = 0; i < N1; ++i){
		a[i] *= b;
	}
}

template<typename T1, std::size_t N1, typename T2>
void operator /=(std::array<T1,N1>& a, const T2& b) noexcept{
	static_assert(std::is_convertible<T2,T1>::value, "Type Mismatch");
	for(std::size_t i = 0; i < N1; ++i){
		a[i] /= b;
	}
}

template<typename T1, std::size_t N1, typename T2>
std::array<T1,N1> operator +(const std::array<T1,N1>& a, const T2& b) noexcept{
	static_assert(std::is_convertible<T2,T1>::value, "Type Mismatch");
	auto ans(a);
	ans += b;
	return ans;
}

template<typename T1, std::size_t N1, typename T2>
std::array<T1,N1> operator +(const T2& b, const std::array<T1,N1>& a) noexcept{
	static_assert(std::is_convertible<T2,T1>::value, "Type Mismatch");
	auto ans(a);
	ans += b;
	return ans;
}

template<typename T1, std::size_t N1, typename T2>
std::array<T1,N1> operator -(const std::array<T1,N1>& a, const T2& b) noexcept{
	static_assert(std::is_convertible<T2,T1>::value, "Type Mismatch");
	auto ans(a);
	ans -= b;
	return ans;
}

template<typename T1, std::size_t N1, typename T2>
std::array<T1,N1> operator -(const T2& b, const std::array<T1,N1>& a) noexcept{
	static_assert(std::is_convertible<T2,T1>::value, "Type Mismatch");
	auto ans(-a);
	ans += b;
	return ans;
}

template<typename T1, std::size_t N1, typename T2>
std::array<T1,N1> operator *(const std::array<T1,N1>& a, const T2& b) noexcept{
	static_assert(std::is_convertible<T2,T1>::value, "Type Mismatch");
	auto ans(a);
	ans *= b;
	return ans;
}

template<typename T1, std::size_t N1, typename T2>
std::array<T1,N1> operator *(const T2& b, const std::array<T1,N1>& a) noexcept{
	static_assert(std::is_convertible<T2,T1>::value, "Type Mismatch");
	auto ans(a);
	ans *= b;
	return ans;
}

template<typename T1, std::size_t N1, typename T2>
std::array<T1,N1> operator /(const std::array<T1,N1>& a, const T2& b) noexcept{
	static_assert(std::is_convertible<T2,T1>::value, "Type Mismatch");
	auto ans(a);
	ans /= b;
	return ans;
}

template<typename T1, std::size_t N1, typename T2>
std::array<T1,N1> operator /(const T2& b, const std::array<T1,N1>& a) noexcept{
	static_assert(std::is_convertible<T2,T1>::value, "Type Mismatch");
	auto ans(a);
	ans /= 1;
	ans *= b;
	return ans;
}


// ============================================================
//                Array / Array Operations
// ============================================================


template<typename T1, std::size_t N1, typename T2, std::size_t N2>
void operator +=(std::array<T1,N1>& a, const std::array<T2,N2>& b) noexcept{
	static_assert(N1 == N2, "Size Mismatch");
	static_assert(std::is_convertible<T2,T1>::value, "Type Mismatch");
	for(std::size_t i = 0; i < N1; ++i){
		a[i] += b[i];
	}
}

template<typename T1, std::size_t N1, typename T2, std::size_t N2>
void operator -=(std::array<T1,N1>& a, const std::array<T2,N2>& b) noexcept{
	static_assert(N1 == N2, "Size Mismatch");
	static_assert(std::is_convertible<T2,T1>::value, "Type Mismatch");
	for(std::size_t i = 0; i < N1; ++i){
		a[i] -= b[i];
	}
}

template<typename T1, std::size_t N1, typename T2, std::size_t N2>
void operator *=(std::array<T1,N1>& a, const std::array<T2,N2>& b) noexcept{
	static_assert(N1 == N2, "Size Mismatch");
	static_assert(std::is_convertible<T2,T1>::value, "Type Mismatch");
	for(std::size_t i = 0; i < N1; ++i){
		a[i] *= b[i];
	}
}

template<typename T1, std::size_t N1, typename T2, std::size_t N2>
void operator /=(std::array<T1,N1>& a, const std::array<T2,N2>& b) noexcept{
	static_assert(N1 == N2, "Size Mismatch");
	static_assert(std::is_convertible<T2,T1>::value, "Type Mismatch");
	for(std::size_t i = 0; i < N1; ++i){
		a[i] /= b[i];
	}
}

template<typename T1, std::size_t N1, typename T2, std::size_t N2>
std::array<T1,N1> operator +(const std::array<T1,N1>& a, const std::array<T2,N2>& b) noexcept{
	static_assert(N1 == N2, "Size Mismatch");
	static_assert(std::is_convertible<T2,T1>::value, "Type Mismatch");
	auto ans(a);
	ans += b;
	return ans;
}

template<typename T1, std::size_t N1, typename T2, std::size_t N2>
std::array<T1,N1> operator -(const std::array<T1,N1>& a, const std::array<T2,N2>& b) noexcept{
	static_assert(N1 == N2, "Size Mismatch");
	static_assert(std::is_convertible<T2,T1>::value, "Type Mismatch");
	auto ans(a);
	ans -= b;
	return ans;
}

template<typename T1, std::size_t N1, typename T2, std::size_t N2>
std::array<T1,N1> operator *(const std::array<T1,N1>& a, const std::array<T2,N2>& b) noexcept{
	static_assert(N1 == N2, "Size Mismatch");
	static_assert(std::is_convertible<T2,T1>::value, "Type Mismatch");
	auto ans(a);
	ans *= b;
	return ans;
}

template<typename T1, std::size_t N1, typename T2, std::size_t N2>
std::array<T1,N1> operator /(const std::array<T1,N1>& a, const std::array<T2,N2>& b) noexcept{
	static_assert(N1 == N2, "Size Mismatch");
	static_assert(std::is_convertible<T2,T1>::value, "Type Mismatch");
	auto ans(a);
	ans /= b;
	return ans;
}








/*
#include <iostream>
template<typename T, std::size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T,N>& a){
	os << "[";
	for(std::size_t i = 0; i < N; ++i){
		if(i) {	os << ","; }
		os << " " << a[i];
	}
	os << "]";
    return os;
}
*/






#endif /* ARRAY_MATH_HPP_ */
