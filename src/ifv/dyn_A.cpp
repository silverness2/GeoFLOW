/*
 * dyn_A.cpp
 *
 *  Created on: Jun 11, 2019
 *      Author: bflynt
 */


#include "dyn_A.hpp"

#include "xstd/array.hpp"


void dyn_A(const IcosGrid& grid, const IcosSoln& soln, IcosSoln& afv){

	using Real    = typename IcosGrid::Real;
	using Integer = typename IcosGrid::Integer;

	const Real g(9.80616);            // Gravity  (m/s^2))
	IcosSoln   Fc(grid.ipe);
	IcosSoln   Fs(grid.ipe);
	//std::vector<std::array<Real,4>> Fc(grid.ipe);
	//std::vector<std::array<Real,4>> Fs(grid.ipe);
	std::vector<Real>               Ec(grid.ipe);
	std::vector<Real>               Es(grid.ipe);

	for(Integer ip = 0; ip < grid.ipe; ++ip){
		const auto nb = grid.nprox[ip];

		const auto vmag2_ip  = dot_product(soln.velo[ip],soln.velo[ip]);

		const auto E0 = 0.5 * vmag2_ip * g * (soln.h[ip] + grid.hb_a[ip]);

	     //s=0.
	     //
	     //-- s1.  h and velocity from hex center to mid edge
	     //
		for(Integer is = 0; is < nb; ++is){
			auto isp = is % nb;
			auto  jp = grid.prox[ip][is];
	        auto  kp = grid.prox[ip][isp];

	        //-- weighted h,vx,vy,vz in Cartesian cord for upper corner

	        Fc.h[is]  = soln.h[ip]*grid.wtsph[ip][is][0]
					  + soln.h[jp]*grid.wtsph[ip][is][1]
					  + soln.h[kp]*grid.wtsph[ip][is][2];
	        for(Integer j = 0; j < IcosSoln::ndim; ++j){
		        Fc.velo[is][j] = soln.velo[ip][j]*grid.wtsph[ip][is][0]
						       + soln.velo[jp][j]*grid.wtsph[ip][is][1]
						       + soln.velo[kp][j]*grid.wtsph[ip][is][2];
	        }
	        //-- interpolate E = 1/2*v*v+g(h+h0) to prepare for grad(E)
	        const auto vmag2_jp = dot_product(soln.velo[jp],soln.velo[jp]);
	        const auto vmag2_kp = dot_product(soln.velo[kp],soln.velo[kp]);

	        Ec[is] = E0 * grid.wtsph[ip][is][0]
				   + 0.5 * vmag2_jp * g * (soln.h[jp] + grid.hb_a[jp]) * grid.wtsph[ip][is][1]
				   + 0.5 * vmag2_kp * g * (soln.h[kp] + grid.hb_a[kp]) * grid.wtsph[ip][is][2];
	    //    s = s + Ec[is];
		}
	     //// E0p=s/dble(nb)  //  E0p = sum ( E_corner + E_corner_m1 ) / 2
	     //// use definition for gradient operator, which is easy to analyze theoretically
		const auto E0p = E0;

	     //--  One-half rule for mid of edge: 1/2 * (X_ims + X_im)
		for(Integer is = 0; is < nb; ++is){
			auto ism = (is+nb-2) % nb;
			Fs.h[is]    = 0.5 * (Fc.h[ism]    + Fc.h[is]);
			Fs.velo[is] = 0.5 * (Fc.velo[ism] + Fc.velo[is]);
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
		for(Integer is = 0; is < nb; ++is){
			auto ism = (is+nb-2) % nb;
			Es[is] = 0.5 * (Ec[ism] + Ec[is]);    // KE on six sides
		}
		Real s, t;
		std::array<Real,3> gradh, gradE;
		for(Integer j = 0; j < 3; ++j){
			s = 0;
			t = 0;
			for(Integer is = 0; is < nb; ++is){
				s += (Fs.h[is] - soln.h[ip]) * grid.Nvec[ip][is][j] * grid.slen[ip][is];
				t += (Es[is] - E0p) * grid.Nvec[ip][is][j] * grid.slen[ip][is];
			}
			gradh[j] = s / grid.area[ip];
			gradE[j] = t / grid.area[ip];
		}

		//
	     //   remove radial component, \vec R component  :  \vec F = \vec F  - (F dot R) \vec R
	     s = dot_product(gradh,grid.Rcvec[ip]);
	     t = dot_product(gradE,grid.Rcvec[ip]);
	     for(Integer j = 0; j < 3; ++j){     // grad = grad - (grad,R) R
	    	 gradh[j] -= s * grid.Rcvec[ip][j];
	    	 gradE[j] -= t * grid.Rcvec[ip][j];
	     }

	     //
	     //-- s3.1 div(V)
	     //-- s3.2 Vorticity = k dot \curl V
	     //
	     Real div_V = 0;
	     Real vort  = 0;
	     for(Integer is = 0; is < nb; ++is){
	    	 s = dot_product(Fs.velo[is], grid.Nvec[ip][is]);   // projection to N
	    	 t = dot_product(Fs.velo[is], grid.Tvec[ip][is]);   // proj to T
	    	 div_V += s * grid.slen[ip][is];
	    	 vort  += t * grid.slen[ip][is];	    	  // (v_i dot T) * ds_i
	     }
	     div_V /= grid.area[ip];
	     vort  /= grid.area[ip];


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
	     for(Integer is = 0; is < nb; ++is){
	    	 s = dot_product(Fs.velo[is], grid.Nvec[ip][is]);   // projection to N
	    	 div_hV += s * Fs.h[is] * grid.slen[ip][is];       // h^*  (V dot N)
	     }
	     div_hV /= grid.area[ip];
	     afv.h[ip] = -div_hV;


	     //
	     //-- s5.  acceleration for velocity
	     //----------------------------------------------------------------------
	     // \p V / \pt + ( xi[eta] + f) * r_unit \cross V = - grad( E = K+gh )
	     //
	     //                          (r_unit \cross V)_i = \Epsil_ijk Rj Vk
	     //----------------------------------------------------------------------
	     //
	     afv.velo[ip] = -(vort+grid.fcori[ip]) * cross_product(grid.Rcvec[ip], soln.velo[ip]) - gradE;

	     //for(Integer i = 0; i < 3; ++i){
	     //	 const auto j = i % 3;
	     //	 const auto k = j % 3;
	     //	 afv.velo[ip][i] = -(vort+grid.fcori[ip]) * (grid.Rcvec[ip][j] * soln.velo[ip][k] - grid.Rcvec[ip][k] * soln.velo[ip][j]) - gradE[i];
	     //}

	     //
	     //-- s6. remove radial component for \part V / \part t
	     //
	     s = dot_product(afv.velo[ip], grid.Rcvec[ip]);
	     afv.velo[ip] -= s * grid.Rcvec[ip];
	}


}
