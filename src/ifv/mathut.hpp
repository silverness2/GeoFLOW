/*
 * mathut.hpp
 *
 *  Created on: Jun 19, 2019
 *      Author: bflynt
 *
 *  This code is a direct translation into C++
 *  from the Fortran code provided by Yonggang Yu.
 */

#ifndef SRC_IFV_MATHUT_HPP_
#define SRC_IFV_MATHUT_HPP_

#include <cmath>

template<typename V, typename M>
void basis_between_sph_car(const V& sph, M& basis){
	auto slon = std::sin(sph[0]);
	auto clon = std::cos(sph[0]);
	auto slat = std::sin(sph[1]);
	auto clat = std::cos(sph[1]);
	basis[0][0] = -slon;               // \vec lambda \dot x  =  - x/rho,  rho=sqrt(x^2+y^2)
	basis[1][0] =  clon;               //             \dot y  =    y/rho
	basis[2][0] =  0;                   //             \dot z  =    0
	basis[0][1] = -slat*clon;      // \vec theta  \dot x  =   -z/r * x/rho
	basis[1][1] = -slat*slon;      //                  y  =   -z/r * y/rho
	basis[2][1] =  clat;               //                  z  =    rho/r
	basis[0][2] =  clat*clon;      // \vec r      \dot x  =    x/r
	basis[1][2] =  clat*slon;      //                  y  =    y/r
	basis[2][2] =  slat;               //                  z  =    z/r
}

template<typename M, typename V>
void AX_mult(const M& A, const V& x, V& b){
	b[0] = A[0][0]*x[0] + A[0][1]*x[1] + A[0][2]*x[2];
	b[1] = A[1][0]*x[0] + A[1][1]*x[1] + A[1][2]*x[2];
	b[2] = A[2][0]*x[0] + A[2][1]*x[1] + A[2][2]*x[2];
}

template<typename V, typename S>
void XY_dot(const V& x, const V& y, S& ans){
	ans = x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}




#endif /* SRC_IFV_MATHUT_HPP_ */
