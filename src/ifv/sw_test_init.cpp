/*
 * sw_test_init.cpp
 *
 *  Created on: Jun 7, 2019
 *      Author: bryan.flynt
 */

#include "sw_test_init.hpp"

#include "mathut.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include "tbox/pio.hpp"

void sw_test_init(const int iswcase, const double alpha, IcosGrid& grid, IcosSoln& soln){
	using Real    = double;
	using Integer = int;
	using std::sin;
	using std::cos;
	using std::acos;
	using std::pow;
	using std::min;
	using std::sqrt;

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

	Real r, ht, h0, c, w, K, R;
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
	//-- s1:  h field defined at center (lat,lon) for both A and C grid
	//
	for(std::size_t ip = 1; ip <= grid.nip; ++ip){
		auto lat = grid.lat(ip);
		auto lon = grid.lon(ip);
		if (iswcase == 1) {
			h0 = 1000;
			r  = acos( sin(lat)*sin(latc) + cos(lat)*cos(latc)*cos(lon-lonc) );
			if (r < R0) {
				ht = 0.5 * h0 * (1 + cos(pi*r/R0));
			}
			else {
				//ht= 0.          // original version
				ht = 0.5 * h0;     // MPAS version
				//ht= 0.1 * h0
			}
		}
		else if (iswcase == 2  || iswcase==3 ) {
			//
			// h0 is 3000 m in case-2
			//
			h0 = 3000;
			// cos\Phi = omega \dot r = cos(delta)
			c  = -sin(alpha)*cos(lat)*cos(lon) + cos(alpha)*sin(lat);
			ht = h0 - (0.5*omegasw*omegasw+omegasw*omega)*ae*ae*c*c/g;
		}
		else if (iswcase == 5) {
			//
			// h0 is 5960 m in case-5
			//
			h0 = 5960;
			// cos\Phi = omega \dot r = cos(delta)
			c  = -sin(alpha)*cos(lat)*cos(lon) + cos(alpha)*sin(lat);
			ht = h0 - (0.5*omegasw*omegasw+omegasw*omega)*ae*ae*c*c/g;
		}
		else if (iswcase == 6) {
			// Rossby-Haurwitz
			h0 = 8000;
			w  = omega_RH;
			K  = K_RH;
			R  = R_RH;
			//
			//        A_RH = 0.5 * w * (2. * omega + w) * cos(lat)**2. + & 0.25 * K**2. * &
			//          cos(lat)**(2.*R) * ((R+1.)*cos(lat)**2. + 2.*R**2. - R - 2. - &
			//          2.*R**2. * cos(lat)**(-2.))
			//
			Real A_RH = 0.5 * w * (2 * omega + w) * cos(lat)*cos(lat)  +
					0.25 * K*K * pow(cos(lat),(2*R)) * ((R+1)*cos(lat)*cos(lat) + 2*R*R - R - 2) -
					0.5  * K*K * pow(cos(lat),(2*R-2)) * R*R;
			//
			Real B_RH = (2*(omega + w)*K / ((R+1)*(R+2))) * pow(cos(lat),R) * ((R*R + 2*R + 2) - ((R+1)*pow(cos(lat),2)));
			//
			Real C_RH = 0.25 * K*K * pow(cos(lat),2*R) * ((R+1.)*cos(lat)*cos(lat) - R - 2);
			//
			ht = h0 + ae*ae/g*(A_RH + B_RH*cos(R*lon) + C_RH*cos(2.*R*lon));
			//        ht = h0 + h0*cos(lat)**2.
		}
		else {
			ht = 0;
		}
		soln(ip,1) = ht;
	}
	//
	//-- terrain surface field
	//
	for(std::size_t ip = 1; ip <= grid.nip; ++ip){
		Real lat = grid.lat(ip);
		Real lon = grid.lon(ip);
		if (iswcase == 5) {
			//
			//  add surface height
			R0 = pi/9.;
			lonc= 1.5*pi;
			latc= pi/6.;
			r  = min(R0, sqrt( pow((lon-lonc),2) + pow((lat-latc),2) ) );
			Real hs = 2000;
			grid.hb_a(ip) = hs * (1. - r/R0);   // surface height
		}
		else {
			grid.hb_a(ip) = 0;
		}
	}
	//
	//note: in C-grid, we will not touch soln(grip.nip+1:niE,1) --- the h-field
	//


	//
	//-- s2:  velocity field
	//
	std::vector<Real> u(grid.NPTS+1);
	std::vector<Real> v(grid.NPTS+1);
	for(std::size_t ipt = 1; ipt <= grid.NPTS; ++ipt){
		Real lat = grid.lat_v(ipt);     // lat velocity
		Real lon = grid.lon_v(ipt);
		Real ut, vt;
		if (iswcase==1 || iswcase==2 || iswcase==5) {
			//
			// geostropic flow
			// V = Omega \cross R
			//
			ut =  u0*(sin(alpha)*sin(lat)*cos(lon) + cos(alpha)*cos(lat));
			vt = -u0*sin(alpha)*sin(lon);
		}
		else if (iswcase == 3) {
			ut =  u0*(cos(lat)*cos(alpha) + sin(lat)*sin(alpha)*cos(lon));
			vt = -u0*sin(alpha)*sin(lon);
		}
		else if (iswcase == 6) {
			ut =  ae*w*cos(lat) + ae*K*pow(cos(lat),R-1)*(R*pow(sin(lat),2) - pow(cos(lat),2))*cos(R*lon);
			vt = -ae*K*R*pow(cos(lat),R-1)*sin(lat)*sin(R*lon);
		}
		else {
			ut = 0;
			vt = 0;
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
	for(std::size_t ipt = 1; ipt <= grid.NPTS; ++ipt){
		Real lat = grid.lat_v(ipt);
		Real lon = grid.lon_v(ipt);
		// unit omega \dot r = cos(delta)
		c = -sin(alpha)*cos(lat)*cos(lon) + cos(alpha)*sin(lat);
		if (iswcase == 1) {
			grid.fcori(ipt)= 0;           // pure advection
		}
		else if (iswcase>=2 || iswcase<=6) {    // case-3 included
			grid.fcori(ipt)= 2.*omega*c;   //  2 * Omega * (\Ome \dot R)
			// global for Earth
		}
		else {
			grid.fcori(ipt)= 0;
		}
	}

	//
	//---  transform velocity (u, v) to cartesian
	//---    (Four vector:  h/Vx/Vy/Vz)
	//
	if (iswcase != 4) {
		std::array<Real,3> P_sph,X,Y;
		std::array<std::array<Real,3>,3> basis;
		for(std::size_t ipt = 1; ipt <= grid.NPTS; ++ipt){
			Real lat = grid.lat_v(ipt);
			Real lon = grid.lon_v(ipt);
			P_sph[0] = lon;
			P_sph[1] = lat;
			P_sph[2] = 1;
			basis_between_sph_car(P_sph,basis);

			X[0] = u[ipt];
			X[1] = v[ipt];
			X[2] = 0;
			AX_mult(basis,X,Y);
			if (grid.stagger == 'A') {
				soln(ipt,2)=Y[0];
				soln(ipt,3)=Y[1];
				soln(ipt,4)=Y[2];
			}
			else if (grid.stagger == 'C') {
				std::array<Real,3> vec1;
				//
				// Proj velocity to N T direction for C-grid
				//      Un=(V,N)N,  Ut=(V,T)T
				//
				auto ip = grid.C4E(1,ipt);
				auto is = grid.C4E(5,ipt);   // Left:  1:4=i,j,k,Lp, 5=is
				//vec1(1:3) = grid.Nvec(1:3,is,ip);     // for iE;  NoutdUn(is,ip)= + 1
				vec1[0] = grid.Nvec(1,is,ip);
				vec1[1] = grid.Nvec(2,is,ip);
				vec1[2] = grid.Nvec(3,is,ip);
				XY_dot(Y,vec1,soln(ipt,2));  //  Un

				vec1[0] = grid.Tvec(1,is,ip);
				vec1[1] = grid.Tvec(2,is,ip);
				vec1[2] = grid.Tvec(3,is,ip);
				XY_dot(Y,vec1,soln(ipt,3));  //  Ut

				//soln(ipt,2) = s;             //  Un
				//soln(ipt,3) = s2;            //  Ut
				//           //
				//           //-- ck
				//           write(6,121) 'iE,is,ip       ', ipt,is,ip
				//           write(6,101) 'grid.Nvec(1:3,is,ip)', grid.Nvec(1:3,is,ip)
				//           write(6,101) 'Tvec(1:3,is,ip)', Tvec(1:3,is,ip)
				//           write(6,101) 'V_xyz(1:3)     ', Y(1:3)
				//           write(6,101) 'VdN            ', s
				//           write(6,101) 'Vdt            ', s2
				//
				// backward test  Un* \vec U_n + Ut * \vec U_T = V = omega \cross R
				// for(std::size_t j = 1; j <= grid.ndim; ++j){
				// 	vec3(j)= s * grid.Nvec(j,is,ip) + s2 * Tvec(j,is,ip);
				// }
			}
			//
			// check with  V = omega \cross R
			//
			//vec1(1) = -sin(alpha);
			//vec1(2) = 0;
			//vec1(3) = cos(alpha);  // \vec Omegasw
			//call vect_from_sph_2_car (P_sph, P_car, ndim);           // \vec R
			//call XY_cross (vec1, P_car, vec2, ndim);                 //  omega \cross R
			//vec2(1:3)=u0*vec2(1:3);
			//call X_2norm(vec2, s2, ndim);
			//call X_2norm(vec3, s3, ndim);
			//call X_2norm(Y, s, ndim);
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
			//        write(6,101) 'Un, Ut    ', soln(ipt,2)  , soln(ipt,3)
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
			for(std::size_t ip = 1; ip <= grid.nip; ++ip){
				ASV_ana(ip, 1)=soln(ip, 1);
			}
			// velocity cartesian
			for(std::size_t ipt = 1; ipt <= grid.NPTS; ++ipt){
				for(std::size_t j = 2; j <= 4; ++j){
					ASV_ana(ipt, j)=soln(ipt, j);
				}
			}
		}
		else if (grid.stagger=='C') {
			for(std::size_t ip = 1; ip <= grid.nip; ++ip){
				ASV_ana(ip, 1)=soln(ip, 1);
			}
			// velocity Un
			for(std::size_t ipt = 1; ipt <= grid.NPTS; ++ipt){
				for(std::size_t j = 2; j <= 3; ++j){
					ASV_ana(ipt, j)=soln(ipt, j);  //  U_n and U_t  vel
				}
			}
		}
	}
	 */

}


