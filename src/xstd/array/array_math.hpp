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






// ============================================================
//                Linear Algebra Operations
// ============================================================


template<typename T1, typename T2, std::size_t N>
typename std::common_type<T1,T2>::type
dot_product(const std::array<T1,N>& a, const std::array<T2,N>& b) noexcept{
	typename std::common_type<T1,T2>::type ans(0);
	for(std::size_t i = 0; i < N; ++i){
		ans += (a[i] * b[i]);
	}
	return ans;
}

template<typename T1, typename T2, std::size_t N>
std::array<typename std::common_type<T1,T2>::type,N>
cross_product(const std::array<T1,N>& a, const std::array<T2,N>& b) noexcept{
	static_assert(3 == N, "Must be Size 3");
	std::array<typename std::common_type<T1,T2>::type,N> ans;
	ans[0] = a[1]*b[2] - a[2]*b[1];
	ans[1] = a[2]*b[0] - a[0]*b[2];
	ans[2] = a[0]*b[1] - a[1]*b[0];
	return ans;
}

template<typename T, std::size_t R, std::size_t C>
std::array<std::array<T,R>,C>
transpose(const std::array<std::array<T,C>,R>& a) noexcept{
	std::array<std::array<T,R>,C> ans;
	for(std::size_t r = 0; r < R; ++r){
		for(std::size_t c = 0; c < C; ++c){
			ans[c][r] = a[r][c];
		}
	}
	return ans;
}

template<typename T1, std::size_t R1, std::size_t C1,
         typename T2, std::size_t R2>
std::array<typename std::common_type<T1,T2>::type,R2>
matmul(const std::array<std::array<T1,C1>,R1>& a,
	   const std::array<T2,R2>& b) noexcept{
	static_assert(C1== R2, "Row/Column Mismatch");
	std::array<typename std::common_type<T1,T2>::type,R2> ans;
	for(std::size_t r = 0; r < R1; ++r){
		ans[r] = 0;
		for(std::size_t i = 0; i < C1; ++i){
			ans[r] += (a[r][i] * b[i]);
		}
	}
	return ans;
}

template<typename T1, std::size_t R1, std::size_t C1,
         typename T2, std::size_t R2, std::size_t C2>
std::array<std::array<typename std::common_type<T1,T2>::type,C2>,R1>
matmul(const std::array<std::array<T1,C1>,R1>& a,
	   const std::array<std::array<T2,C2>,R2>& b) noexcept{
	static_assert(C1== R2, "Row/Column Mismatch");
	std::array<std::array<typename std::common_type<T1,T2>::type,C2>,R1> ans;
	for(std::size_t r = 0; r < R1; ++r){
		for(std::size_t c = 0; c < C2; ++c){
			ans[r][c] = 0;
			for(std::size_t i = 0; i < C1; ++i){
				ans[r][c] += (a[r][i] * b[i][c]);
			}
		}
	}
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
