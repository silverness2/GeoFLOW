/*
 * dyn_A.cpp
 *
 *  Created on: Jun 11, 2019
 *      Author: bflynt
 *
 *  NOTE:
 *	This code is a C++ representation of Fortran code
 *	provided by Yonggang Yu (CIRES/NOAA) to solve the
 *	2-D Shallow Water Equations on a sphere.
 */


#include "dyn_A.hpp"

#include <array>
#include <cmath>
#include <vector>


void dyn_A(const int iswcase, const IcosGrid& grid, const IcosSoln& soln, IcosSoln& afv){

	using Real    = typename IcosGrid::Real;
	using Integer = typename IcosGrid::Integer;

	const Real g(9.80616);            // Gravity  (m/s^2))
	IcosSoln   Fc,Fs;
	Fc.resize(4,grid.NPTS).reorder(IcosSoln::f_order).rebase(1);
	Fs.resize(4,grid.NPTS).reorder(IcosSoln::f_order).rebase(1);
	std::vector<Real> Ec(grid.ipe+1);
	std::vector<Real> Es(grid.ipe+1);

	for(Integer ip = grid.ips; ip <= grid.ipe; ++ip){
		const auto nb = grid.nprox(ip);

		const auto vmag2_ip  = std::pow(soln(ip,2),2) + std::pow(soln(ip,3),2) + std::pow(soln(ip,4),2);
		const auto E0 = 0.5 * vmag2_ip * g * (soln(ip,1) + grid.hb_a(ip));

	     //s=0.
	     //
	     //-- s1.  h and velocity from hex center to mid edge
	     //
		for(Integer is = 1; is <= nb; ++is){
			auto isp = (is % nb) + 1;
			auto  jp = grid.prox(is,ip);
	        auto  kp = grid.prox(isp,ip);

	        //-- weighted h,vx,vy,vz in Cartesian cord for upper corner
	        for(Integer j = 1; j <= grid.ndimFV; ++j){
		        Fc(j,is) = soln(ip,j)*grid.wtsph(1,is,ip)
						 + soln(jp,j)*grid.wtsph(2,is,ip)
						 + soln(kp,j)*grid.wtsph(3,is,ip);
	        }
	        //-- interpolate E = 1/2*v*v+g(h+h0) to prepare for grad(E)
	        const auto vmag2_jp  = std::pow(soln(jp,2),2) + std::pow(soln(jp,3),2) + std::pow(soln(jp,4),2);
	        const auto vmag2_kp  = std::pow(soln(kp,2),2) + std::pow(soln(kp,3),2) + std::pow(soln(kp,4),2);


	        Ec[is] = E0 * grid.wtsph(1,is,ip)
				   + 0.5 * vmag2_jp * g * (soln(jp,1) + grid.hb_a(jp)) * grid.wtsph(2,is,ip)
				   + 0.5 * vmag2_kp * g * (soln(kp,1) + grid.hb_a(kp)) * grid.wtsph(3,is,ip);
	    //    s = s + Ec[is];
		}
	     //// E0p=s/dble(nb)  //  E0p = sum ( E_corner + E_corner_m1 ) / 2
	     //// use definition for gradient operator, which is easy to analyze theoretically
		const auto E0p = E0;

	     //--  One-half rule for mid of edge: 1/2 * (X_ims + X_im)
		for(Integer is = 1; is <= nb; ++is){
			auto ism = ((is+nb-2) % nb) + 1;
			for(Integer j = 1; j <= grid.ndimFV; ++j){
				Fs(j,is)= 0.5*(Fc(j,ism) + Fc(j,is));
			}
		}

	     //
	     //-- s2.1 gradient h
	     //-- s2.2 Kinetic energy Gradient :  K + g (h + hb_a)
	     //        \grad E_j = 1/A  *  \sum_i (Es_i -E0) * n(j,i)  dsi
	     //-- s2.3 remove radial component
	     //
	     //
	     //-- use soln.h[ip]
	     ////havg=0.
	     ////do is=1,nb
	     ////   havg= havg + Fs(1,is)
	     ////enddo
	     ////havg= havg/dble(nb)
	     //
		for(Integer is = 1; is <= nb; ++is){
			auto ism = ((is+nb-2) % nb) + 1;
			Es[is] = 0.5 * (Ec[ism] + Ec[is]);    // KE on six sides
		}
		Real s, t;
		std::array<Real,3> gradh, gradE;
		for(Integer j = 1; j <= 3; ++j){
			s = 0;
			t = 0;
			for(Integer is = 1; is <= nb; ++is){
				s += (Fs(1,is) - soln(ip,1)) * grid.Nvec(j,is,ip) * grid.slen(is,ip); // (h-hc)* n ds; 1=h, 2_4=Vx_z
				t += (Es[is] - E0p)*grid.Nvec(j,is,ip)*grid.slen(is,ip);
			}
			gradh[j-1] = s / grid.area(ip);
			gradE[j-1] = t / grid.area(ip);
		}

		//
	    //   remove radial component, \vec R component  :  \vec F = \vec F  - (F dot R) \vec R
		s = gradh[0]*grid.Rcvec(1,ip) + gradh[1]*grid.Rcvec(2,ip) + gradh[2]*grid.Rcvec(3,ip);
		t = gradE[0]*grid.Rcvec(1,ip) + gradE[1]*grid.Rcvec(2,ip) + gradE[2]*grid.Rcvec(3,ip);
	     for(Integer j = 1; j <= 3; ++j){     // grad = grad - (grad,R) R
	    	 gradh[j-1] -= s * grid.Rcvec(j,ip);
	    	 gradE[j-1] -= t * grid.Rcvec(j,ip);
	     }

	     //
	     //-- s3.1 div(V)
	     //-- s3.2 Vorticity = k dot \curl V
	     //
	     Real div_V = 0;
	     Real vort  = 0;
	     for(Integer is = 1; is <= nb; ++is){
	    	 s = 0;
	    	 t = 0;
	    	 for(Integer j = 1; j <= 3; ++j){
	    		 s += Fs(j+1,is)*grid.Nvec(j,is,ip);
	    		 t += Fs(j+1,is)*grid.Tvec(j,is,ip);
	    	 }
	    	 div_V += s * grid.slen(is,ip);
	    	 vort  += t * grid.slen(is,ip);	    	  // (v_i dot T) * ds_i
	     }
	     div_V /= grid.area(ip);
	     vort  /= grid.area(ip);


	//     //
	//     //-- s4.  div (hV) =  < V(center), \grad (h) >  + h(center)* div_V
	//     //
	//     s=0.               // s = (V, g(h))
	//     do j=1,3
	//        s= s + gradh(j)*FV(ip,j+1)
	//     enddo
	//     afv(ip,1)= -(FV(ip,1)*div_V + s)


	     //
	     //-- s4.  div (hV) =  \int  (hV \dot N) * ds  / A
	     //
	     Real div_hV = 0;
	     for(Integer is = 1; is <= nb; ++is){
	    	 s = 0;
	    	 for(Integer j = 1; j <= 3; ++j){
	    		 s += Fs(j+1,is)*grid.Nvec(j,is,ip); // projection to N
	    	 }
	    	 div_hV += s*Fs(1,is)*grid.slen(is,ip);  // h^*  (V dot N)
	     }
	     div_hV /= grid.area(ip);
	     afv(ip,1) = -div_hV;


	     //
	     //-- s5.  acceleration for velocity
	     //----------------------------------------------------------------------
	     // \p V / \pt + ( xi[eta] + f) * r_unit \cross V = - grad( E = K+gh )
	     //
	     //                          (r_unit \cross V)_i = \Epsil_ijk Rj Vk
	     //----------------------------------------------------------------------
	     //
	     for(Integer i = 1; i <= 3; ++i){
	    	 auto j = (i % 3) + 1;
	    	 auto k = (j % 3) + 1;
	    	 afv(ip,i+1) = -(vort+grid.fcori(ip))*(grid.Rcvec(j,ip)*soln(ip,k+1)-grid.Rcvec(k,ip)*soln(ip,j+1))-gradE[i-1];
	     }

	     //for(Integer i = 0; i < 3; ++i){
	     //	 const auto j = i % 3;
	     //	 const auto k = j % 3;
	     //	 afv.velo[ip][i] = -(vort+grid.fcori[ip]) * (grid.Rcvec[ip][j] * soln.velo[ip][k] - grid.Rcvec[ip][k] * soln.velo[ip][j]) - gradE[i];
	     //}

	     //
	     //-- s6. remove radial component for \part V / \part t
	     //
	     s = 0;
	     for(Integer j = 1; j <= 3; ++j){
	    	 s += afv(ip,j+1) * grid.Rcvec(j,ip);
	     }
	     for(Integer j = 1; j <= 3; ++j){
	    	 afv(ip,j+1) = afv(ip,j+1) - s*grid.Rcvec(j,ip);
	     }
	}


	  //
	  //------------------------------------------------------------
	  //   reset acceleration if sw_test_case_1
	  //         passive advection only
	  //------------------------------------------------------------
	  //
	  if ( iswcase == 1 ) {
			for(Integer ip = grid.ips; ip <= grid.ipe; ++ip){
				for(Integer j = 2; j <= 4; ++j){
					afv(ip,j) = 0;
				}
			}
	  }


}
