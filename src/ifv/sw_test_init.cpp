/*
 * sw_test_init.cpp
 *
 *  Created on: Jun 7, 2019
 *      Author: bryan.flynt
 */

#include "sw_test_init.hpp"

#include "xstd/array.hpp"

#include <algorithm>
#include <array>
#include <cmath>


void sw_test_init(const int iswcase, const double alpha, IcosGrid& grid, IcosSoln& soln){
	using Real    = double;
	using Integer = int;

	//
	// dimension: meter, m/s
	//
	const Real pi(3.1415926535897932384626433);
	const Real ae(6371220.0);         // Earth Radius (m)
	const Real g(9.80616);            // Gravity  (m/s^2))
	const Real omega(7.292e-5);       // earth rotation rate (rad/s)
	const Real omega_RH(7.848e-6);    // s^-1
	const Real K_RH(7.848e-6);        // s^-1
	const Real R_RH(4.0);
	const Real thetab(-pi/6.0);       // TC3
	const Real thetae(pi/2.0);
	const Real xe(0.30);

	Real ht,c;
	Real R0      = 1;
	Real latc    = 0;
	Real lonc    = -pi/2.0;
	Real omegasw = 2.0*pi/(12.0*86400.0); // 2pi/T
	Real u0      = ae*omegasw;
	if (iswcase == 4) {
		u0 =  20.0;             // 20 or 40 m/s
		omegasw= u0/ae;         // redefine
	}
	else if (iswcase==5) {
		u0=  20.0;            // 20m/s
	}


	//
	//-- s1:  h field defined at center (lat_h,lon_h) for both A and C grid
	//
	for(Integer ip = 0; ip < grid.nip; ++ip) {
		const auto lat= grid.lat[ip];
		const auto lon= grid.lon[ip];
		if (iswcase == 1) {
			const Real h0 = 1000.0;
			const Real r = std::acos( std::sin(lat)*std::sin(latc) + std::cos(lat)*std::cos(latc)*std::cos(lon-lonc) );
			if (r < R0) {
				ht = 0.50 * h0 * (1.0 + std::cos(pi*r/R0));
			}
			else {
				//ht= 0.0          // original version
				ht = 0.50 * h0;     // MPAS version
				//ht= 0.10 * h0
			}
		}
		else if (iswcase == 2  || iswcase==3 ) {
			//
			// h0 is 3000 m in case-2
			//
			const Real h0= 3000.0;
			// cos\Phi = omega \dot r = std::cos(delta);
			c = -std::sin(alpha)*std::cos(lat)*std::cos(lon) + std::cos(alpha)*std::sin(lat);
			ht = h0 - (0.50*omegasw*omegasw+omegasw*omega)*ae*ae*c*c/g;
		}
		else if (iswcase == 5) {
			//
			// h0 is 5960 m in case-5
			//
			const Real h0= 5960.0;
			// cos\Phi = omega \dot r = std::cos(delta)
			c  = -std::sin(alpha)*std::cos(lat)*std::cos(lon) + std::cos(alpha)*std::sin(lat);
			ht = h0 - (0.50*omegasw*omegasw+omegasw*omega)*ae*ae*c*c/g;
		}
		else if (iswcase == 6) {
			// Rossby-Haurwitz
			const Real h0 = 8.0e+3;
			const Real w = omega_RH;
			const Real K = K_RH;
			const Real R = R_RH;
			//
			//        A_RH = 0.50 * w * (2.0 * omega + w) * std::cos(lat)**2.0 + & 0.250 * K**2.0 * &
			//          std::cos(lat)**(2.0*R) * ((R+1.0)*std::cos(lat)**2.0 + 2.0*R**2.0 - R - 2.0 - &
			//          2.0*R**2.0 * std::cos(lat)**(-2.0))
			//
			const Real A_RH =
					0.50 * w * (2.0 * omega + w) * std::cos(lat)*std::cos(lat) +
					0.250 * K*K * std::pow(std::cos(lat),(2.0*R)) * ((R+1.0)*std::cos(lat)*std::cos(lat) + 2.0*R*R - R - 2.0) -
					0.50  * K*K * std::pow(std::cos(lat),(2.0*R-2.0)) * R*R;
			//
			const Real B_RH =
					(2.0*(omega + w)*K / ((R+1.0)*(R+2.0))) * std::pow(std::cos(lat),R) * ((R*R +
							2.0*R + 2.0) - std::pow(((R+1.0)*std::cos(lat)),2.0));
			//
			const Real C_RH =
					0.250 * K*K * std::pow(std::cos(lat),(2.0*R)) * ((R+1.0)*std::cos(lat)*std::cos(lat) - R - 2.0);
			//
			ht = h0 + ae*ae/g*(A_RH + B_RH*std::cos(R*lon) + C_RH*std::cos(2.0*R*lon));
			//        ht = h0 + h0*std::cos(lat)**2.0;
		}
		else {
			ht = 0.0;
		}
		soln.h[ip] = ht;
	}

	//
	//-- terrain surface field
	//
	for(Integer ip = 0; ip < grid.nip; ++ip){
		const auto lat= grid.lat[ip];
		const auto lon= grid.lon[ip];
		if (iswcase == 5) {
			//
			//  add surface height
			R0 = pi/9.0;
			lonc= 1.50*pi;
			latc= pi/6.0;
			const Real r = std::min(R0, std::sqrt( std::pow((lon-lonc),2) + std::pow(lat-latc,2)));
			const Real hs = 2000.0;
			grid.hb_a[ip]= hs * (1.0 - r/R0);  // surface height
		}
		else {
			grid.hb_a[ip]= 0.0;
		}
	}
	//
	//note: in C-grid, we will not touch FV(nip+1:niE,1) --- the h-field
	//


	//
	//-- s2:  velocity field
	//
	Real ut, vt;
	std::vector<Real> u(grid.NPTS);
	std::vector<Real> v(grid.NPTS);
	for(Integer ipt = 0; ipt < grid.NPTS; ++ipt){
		const auto lat = grid.lat_v[ipt];     // lat velocity
		const auto lon = grid.lon_v[ipt];     //
		if (iswcase==1 || iswcase==2 || iswcase==5) {
			//
			// geostropic flow
			// V = Omega \cross R
			//
			ut =  u0*(std::sin(alpha)*std::sin(lat)*std::cos(lon) + std::cos(alpha)*std::cos(lat));
			vt = -u0*std::sin(alpha)*std::sin(lon);
		}
		else if (iswcase == 3) {
			ut =  u0*(std::cos(lat)*std::cos(alpha) + std::sin(lat)*std::sin(alpha)*std::cos(lon));
	        vt = -u0*std::sin(alpha)*std::sin(lon);
		}
		else if (iswcase == 6) {
			const Real w = omega_RH;
			const Real K = K_RH;
			const Real R = R_RH;
			ut =  ae*w*std::cos(lat) + std::pow(ae*K*std::cos(lat),(R-1.0))*(   std::pow(R*std::sin(lat),2.0) - std::pow(std::cos(lat),2.0))*std::cos(R*lon);
			vt = -ae*K*R*std::pow(std::cos(lat),(R-1.0))*std::sin(lat)*std::sin(R*lon);
		}
		else {
			ut = 0.0;
			vt = 0.0;
		}
		u[ipt] = ut;
		v[ipt] = vt;
	}


	//
	//-- s3:
	//--      Coriolis (general for both A and C), unified
	//        note velocity and acceleration
	//        for A-grid at center, C-grid on edge
	//
	for(Integer ipt = 0; ipt < grid.NPTS; ++ipt){
		const auto lat = grid.lat_v[ipt];
		const auto lon = grid.lon_v[ipt];
		// unit omega \dot r = std::cos(delta)
		c  = -std::sin(alpha)*std::cos(lat)*std::cos(lon) + std::cos(alpha)*std::sin(lat);
		if (iswcase == 1) {
			grid.fcori[ipt]= 0.0;          // pure advection
		}
		else if (iswcase >= 2 || iswcase <= 6) {    // case-3 included
			grid.fcori[ipt]= 2.0*omega*c;   //  2 * Omega * (\Ome \dot R)
		}                              // global for Earth
		else {
			grid.fcori[ipt]= 0.0;
		}
	}


	//
	//---  transform velocity (u, v) to cartesian
	//---    (Four vector:  h/Vx/Vy/Vz)
	//
	if (iswcase != 4) {
		for(Integer ipt = 0; ipt < grid.NPTS; ++ipt){
			const auto lat = grid.lat_v[ipt];
			const auto lon = grid.lon_v[ipt];

			const auto clat = std::cos(lat);
			const auto slat = std::sin(lat);
			const auto clon = std::cos(lon);
			const auto slon = std::sin(lon);
			std::array<Real,IcosGrid::ndim> P_sph({lon,lat,0});
			std::array<std::array<Real,IcosGrid::ndim>,IcosGrid::ndim> basis;
			basis[0] = { -slon, -slat*clon, clat*clon};
			basis[1] = {  clon, -slat*slon, clat*slon};
			basis[2] = {     0,       clat,      slat};

			std::array<Real,IcosGrid::ndim> X({u[ipt],v[ipt],0});

			auto Y = matmul(basis,X);    // Y(1:3) cartesian velocity
			if (grid.stagger == 'A') {
				soln.velo[ipt] = Y;
				//FV.velo[ipt][0] = Y[0];
				//FV.velo[ipt][1] = Y[1];
				//FV.velo[ipt][2] = Y[2];
				//vec3(1:3)=0.0;       //  aux
			}
			else if (grid.stagger == 'C') {
				//
				// Proj velocity to N T direction for C-grid
				//      Un=(V,N)N,  Ut=(V,T)T
				//
				const auto ip   = grid.C4E[ipt][0]; //  (1,ipt);
				const auto is   = grid.C4E[ipt][4]; //  (5,ipt); // Left:  1:4=i,j,k,Lp, 5=is
				const auto vec1 = grid.Nvec[ip][is];
				auto s = dot_product(Y,vec1);

				const auto vec2 = grid.Tvec[ip][is];
				auto s2 = dot_product(Y,vec2);

				soln.velo[ipt][0]=s;             //  Un
				soln.velo[ipt][1]=s2;            //  Ut
				soln.velo[ipt][2]=0;
				//           //
				//           //-- ck
				//           write(6,121) 'iE,is,ip       ', ipt,is,ip
				//           write(6,101) 'Nvec(1:3,is,ip)', Nvec(1:3,is,ip)
				//           write(6,101) 'Tvec(1:3,is,ip)', Tvec(1:3,is,ip)
				//           write(6,101) 'V_xyz(1:3)     ', Y(1:3)
				//           write(6,101) 'VdN            ', s
				//           write(6,101) 'Vdt            ', s2
				//
				// backward test  Un* \vec U_n + Ut * \vec U_T = V = omega \cross R
				// for(Integer j = 0; j < 3; ++j){
				// 	vec3(j)= s * Nvec(j,is,ip) + s2 * Tvec(j,is,ip);
				// }
			}
			//
			// check with  V = omega \cross R
			//
			// vec1(1)=-std::sin(alpha);
			// vec1(2)=0.0;
			// vec1(3)=std::cos(alpha);  // \vec Omegasw
			// call vect_from_sph_2_car (P_sph, P_car, IcosGrid::ndim);           // \vec R
			// call XY_cross (vec1, P_car, vec2, IcosGrid::ndim);                 //  omega \cross R
			// vec2(1:3)=u0*vec2(1:3);
			// call X_2norm(vec2, s2, IcosGrid::ndim);
			// call X_2norm(vec3, s3, IcosGrid::ndim);
			// call X_2norm(Y, s, IcosGrid::ndim);
			//        //
			//        //-- ck
			//        //
			//        write(6,121) 'ipt=      ', ipt
			//        write(6,101) 'lon, lat  ', lon, lat
			//        write(6,101) '\vec Omega', vec1(1:3)
			//        write(6,101) '\vec R    ', P_car(1:3)
			//        write(6,101) '\vec Omega \cross R', vec2(1:3)
			//        write(6,101) 'vec1(1:3)-Y(1:3), norm diff', vec2(1:3)-Y(1:3), s-s2
			//        write(6,101) 'vec3(1:3)-Y(1:3), norm diff', vec3(1:3)-Y(1:3), s-s3
			//        write(6,101) 'Un, Ut    ', FV(ipt,2)  , FV(ipt,3)
			//        write(6,'(//)')
			//
		}      // ipt
	}


/*
	if (iswcase==2 || iswcase==3) {
		// steady state let's save initial state to compute Error norms
		// 1  2  3  4
		// h  Vx Vy Vz
		//
		// h
		if (grid.stagger=='A') {
			for(Integer ip = 0; ip < grid.nip; ++ip){
				ASV_ana(ip, 0) = FV[ip][0];
			}
			// velocity cartesian
			for(Integer ipt = 0; ipt < NPTS; ++ipt){
				for(Integer j = 1; j < 4; ++j){
					ASV_ana(ipt, j) = FV[ipt][j];
				}
			}
		}
		else if (grid.stagger=='C') {
			for(Integer ip = 0; ip < grid.nip; ++ip){
				ASV_ana(ip, 0) = FV[ip][0];
			}
			// velocity Un
			for(Integer ipt = 0; ipt < NPTS; ++ipt){
				for(Integer j = 1; j < 3; ++j){
					ASV_ana(ipt, j) = FV[ipt][j];  //  U_n and U_t  vel
				}
			}
		}
	}
*/



}


