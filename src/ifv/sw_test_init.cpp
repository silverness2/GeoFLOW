/*
 * sw_test_init.cpp
 *
 *  Created on: Jun 7, 2019
 *      Author: bryan.flynt
 */

#include "sw_test_init.hpp"

#include <cmath>

void sw_test_init(IcosGrid& grid){
	using Real    = double;
	using Integer = int;

	//
	// dimension: meter, m/s
	//
	const Real pi(3.1415926535897932384626433);
	Real R0      = 1;
	Real latc    = 0;
	Real lonc    = -pi/2.0;
	Real omegasw = 2.0*pi/(12.0*86400.0); // 2pi/T
	if (iswcase == 4) {
		u0 =  20.0;             // 20 or 40 m/s
		omegasw= u0/ae         // redefine
	}
	else if (iswcase==5) {
		u0=  20.0;            // 20m/s
	}
	else {
		u0=  ae*omegasw;
	}
	ndim=    3;                        // fix
	nx=      nx_NICAM;
	mx=      mx_NICAM;
	omega_RH=7.848e-6;                 // s^-1
	K_RH    =7.848e-6;                 // s^-1
	R_RH    =4.0;                     //
	thetab=-pi/6.0;                   // TC3
	thetae= pi/2.0;
	xe= 0.30;

	//
	//-- s1:  h field defined at center (lat_h,lon_h) for both A and C grid
	//
	for(Integer ip = 0; ip < nip; ++ip) {
		lat= lat_h(ip);
		lon= lon_h(ip);
		if (iswcase == 1) {
			h0=1000.0;
			r  = std::acos( std::sin(lat)*std::sin(latc) + std::cos(lat)*std::cos(latc)*std::cos(lon-lonc) );
			if (r < R0) {
				ht= 0.50 * h0 * (1.0 + std::cos(pi*r/R0));
			}
			else {
				//ht= 0.0          // original version
				ht= 0.50 * h0;     // MPAS version
				//ht= 0.10 * h0
			}
		}
		else if (iswcase == 2  || iswcase==3 ) {
			//
			// h0 is 3000 m in case-2
			//
			h0= 3000.0;
			// cos\Phi = omega \dot r = std::cos(delta);
			c  = -std::sin(alpha)*std::cos(lat)*std::cos(lon) + std::cos(alpha)*std::sin(lat);
			ht = h0 - (0.50*omegasw*omegasw+omegasw*omega)*ae*ae*c*c/g;
		}
		else if (iswcase == 5) {
			//
			// h0 is 5960 m in case-5
			//
			h0= 5960.0;
			// cos\Phi = omega \dot r = std::cos(delta)
			c  = -std::sin(alpha)*std::cos(lat)*std::cos(lon) + std::cos(alpha)*std::sin(lat);
			ht = h0 - (0.50*omegasw*omegasw+omegasw*omega)*ae*ae*c*c/g;
		}
		else if (iswcase == 6) {
			// Rossby-Haurwitz
			h0= 8.d3;
			w = omega_RH;
			k = K_RH;
			R = R_RH;
			//
			//        A_RH = 0.50 * w * (2.0 * omega + w) * std::cos(lat)**2.0 + & 0.250 * K**2.0 * &
			//          std::cos(lat)**(2.0*R) * ((R+1.0)*std::cos(lat)**2.0 + 2.0*R**2.0 - R - 2.0 - &
			//          2.0*R**2.0 * std::cos(lat)**(-2.0))
			//
			A_RH = 0.50 * w * (2.0 * omega + w) * std::cos(lat)*std::cos(lat)  +
					0.250 * K*K * std::pow(std::cos(lat),(2.0*R)) * ((R+1.0)*std::cos(lat)*std::cos(lat) + 2.0*R*R - R - 2.0) -
					0.50  * K*K * std::pow(std::cos(lat),(2.0*R-2.0)) * R*R;
			//
			B_RH = (2.0*(omega + w)*K / ((R+1.0)*(R+2.0))) * std::pow(std::cos(lat),R) * ((R*R +
					2.0*R + 2.0) - std::pow(((R+1.0)*std::cos(lat)),2.0));
			//
			C_RH = 0.250 * K*K * std::pow(std::cos(lat),(2.0*R)) * ((R+1.0)*std::cos(lat)*std::cos(lat) - R - 2.0);
			//
			ht = h0 + ae*ae/g*(A_RH + B_RH*std::cos(R*lon) + C_RH*std::cos(2.0*R*lon));
			//        ht = h0 + h0*std::cos(lat)**2.0;
		}
		else
			ht=0.0;
	}
	FV(ip,1)=ht;

	//
	//-- terrain surface field
	//
	for(Integer ip = 0; ip < nip; ++ip){
		lat= lat_h(ip);
		lon= lon_h(ip);
		if (iswcase == 5) {
			//
			//  add surface height
			R0 = pi/9.0;
			lonc= 1.50*pi;
			latc= pi/6.0;
			r  = min(R0, std::sqrt( std::pow((lon-lonc),2) + std::pow(lat-latc,2)));
			hs = 2000.0;
			hb_a(ip)= hs * (1.0 - r/R0);  // surface height
		}
		else {
			hb_a(ip)= 0.0;
		}
	}
	//
	//note: in C-grid, we will not touch FV(nip+1:niE,1) --- the h-field
	//


	//
	//-- s2:  velocity field
	//
	for(Integer ipt = 0; ipt < NPTS; ++ipt){
		lat= lat_v(ipt);     // lat velocity
		lon= lon_v(ipt);     //
		if (iswcase==1 || iswcase==2 || iswcase==5) {
			//
			// geostropic flow
			// V = Omega \cross R
			//
			ut =  u0*(std::sin(alpha)*std::sin(lat)*std::cos(lon) + std::cos(alpha)*std::cos(lat));
			vt = -u0*std::sin(alpha)*std::sin(lon);
		}
		else if (iswcase == 3) {
			ut =  u0*(std::cos(lat)*std::cos(alpha) + std::sin(lat)*std::sin(alpha)*std::cos(lon))
	        						vt = -u0*std::sin(alpha)*std::sin(lon);
		}
		else if (iswcase == 6) {
			ut =  ae*w*std::cos(lat) + ae*K*std::cos(lat)**(R-1.0)*(R*std::sin(lat)**2.0 - std::cos(lat)**2.0)*std::cos(R*lon);
			vt = -ae*K*R*std::cos(lat)**(R-1.0)*std::sin(lat)*std::sin(R*lon);
		}
		else {
			ut = 0.0;
			vt = 0.0;
		}
		u(ipt)=ut; v(ipt)=vt;
	}


	//
	//-- s3:
	//--      Coriolis (general for both A and C), unified
	//        note velocity and acceleration
	//        for A-grid at center, C-grid on edge
	//
	for(Integer ipt = 0; ipt < NPTS; ++ipt){
		lat= lat_v(ipt);
		lon= lon_v(ipt);
		// unit omega \dot r = std::cos(delta)
		c  = -std::sin(alpha)*std::cos(lat)*std::cos(lon) + std::cos(alpha)*std::sin(lat);
		if (iswcase == 1) {
			fcori(ipt)= 0.0;          // pure advection
		}
		else if (iswcase.GE.2 || iswcase.LE.6) {    // case-3 included
			fcori(ipt)= 2.0*omega*c;   //  2 * Omega * (\Ome \dot R)
		}                              // global for Earth
		else {
			fcori(ipt)= 0.0;
		}
	}


	//
	//---  transform velocity (u, v) to cartesian
	//---    (Four vector:  h/Vx/Vy/Vz)
	//
	if (iswcase != 4) {
		for(Integer ipt = 0; ipt < NPTS; ++ipt){
			lat= lat_v(ipt);
			lon= lon_v(ipt);
			P_sph(1)=lon;
			P_sph(2)=lat;
			P_sph(3)=1.0;
			call basis_between_sph_car (P_sph, basis, invbasis, ndim)
			X(1)=u(ipt);
			X(2)=v(ipt);
			X(3)=0.0;
			m=ndim;
			n=ndim;
			mp=m;
			np=n;
			call AX_mult(basis, X, Y, m, n, mp, np)       // Y(1:3) cartesian velocity
					if (stagger == 'A') {
						FV(ipt,2)=Y(1);
						FV(ipt,3)=Y(2);
						FV(ipt,4)=Y(3);
						vec3(1:3)=0.0;       //  aux
					}
					else if (stagger == 'C') {
						//
						// Proj velocity to N T direction for C-grid
						//      Un=(V,N)N,  Ut=(V,T)T
						//
						ip=C4E(1,ipt);
						is=C4E(5,ipt);   // Left:  1:4=i,j,k,Lp, 5=is
						vec1(1:3)= Nvec(1:3,is,ip);     // for iE;  NoutdUn(is,ip)= + 1
						call XY_dot (Y,vec1,n,ndim,s);
						vec1(1:3)=  Tvec(1:3,is,ip);   //
						call XY_dot (Y,vec1,n,ndim,s2);
						FV(ipt,2)=s;             //  Un
						FV(ipt,3)=s2;            //  Ut
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
						for(Integer j = 0; j < 3; ++j){
							vec3(j)= s * Nvec(j,is,ip) + s2 * Tvec(j,is,ip);
						}
					}
			//
			// check with  V = omega \cross R
			//
			vec1(1)=-std::sin(alpha);
			vec1(2)=0.0;
			vec1(3)=std::cos(alpha);  // \vec Omegasw
			call vect_from_sph_2_car (P_sph, P_car, ndim);           // \vec R
			call XY_cross (vec1, P_car, vec2, ndim);                 //  omega \cross R
			vec2(1:3)=u0*vec2(1:3);
			call X_2norm(vec2, s2, ndim);
			call X_2norm(vec3, s3, ndim);
			call X_2norm(Y, s, ndim);
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



	if (iswcase==2 || iswcase==3) {
		// steady state let's save initial state to compute Error norms
		// 1  2  3  4
		// h  Vx Vy Vz
		//
		// h
		if (stagger=='A') {
			for(Integer ip = 0; ip < nip; ++ip){
				ASV_ana(ip, 1)=FV(ip, 1);
			}
			// velocity cartesian
			for(Integer ipt = 0; ipt < NPTS; ++ipt){
				for(Integer i = 1; i < 4; ++i){
					ASV_ana(ipt, j)=FV(ipt, j);
				}
			}
		}
		else if (stagger=='C') {
			for(Integer ip = 0; ip < nip; ++ip){
				ASV_ana(ip, 1)=FV(ip, 1);
			}
			// velocity Un
			for(Integer ipt = 0; ipt < NPTS; ++ipt){
				for(Integer i = 1; i < 3; ++i){
					ASV_ana(ipt, j)=FV(ipt, j);  //  U_n and U_t  vel
				}
			}
		}
	}



}


