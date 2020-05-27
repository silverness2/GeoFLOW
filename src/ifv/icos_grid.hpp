/*
 * grid.hpp
 *
 *  Created on: Jun 6, 2019
 *      Author: bflynt
 *
 *  This code is a direct translation into C++
 *  from the Fortran code provided by Yonggang Yu.
 */

#ifndef SRC_IFV_GRID_HPP_
#define SRC_IFV_GRID_HPP_

//#include <array>
#include <cstdint>
#include <string>
//#include <vector>

#include "tbox/multi_array.hpp"


struct IcosGrid {

public:
	using Real    = double;
	using Integer = std::int32_t;

	IcosGrid();
	IcosGrid(const IcosGrid& other) = default;
	IcosGrid(IcosGrid&& other) = default;
	~IcosGrid() = default;

	IcosGrid& operator=(const IcosGrid& other) = default;
	IcosGrid& operator=(IcosGrid&& other) = default;

	void load(const char stagger, const Integer level, const std::string filename);


public:
	static constexpr Integer ndim  = 3;
	static constexpr Integer nside = 6;

	/*
	std::vector<std::array<Real,2>>                      icos_grid;
	std::vector<std::array<std::array<Real,ndim>,nside>> wtsph;
	std::vector<std::array<Real,nside>>                  slen;
	std::vector<Real>                                    area;
	std::vector<std::array<std::array<Real,ndim>,nside>> Nvec;
	std::vector<std::array<std::array<Real,ndim>,nside>> Tvec;
	std::vector<std::array<std::array<Real,ndim>,nside>> Rvec;
	std::vector<std::array<Real,ndim>>                   Rcvec;
	std::vector<Integer>                                 nprox;
	std::vector<std::array<Integer,nside>>               prox;
	std::vector<std::array<Integer,nside>>               proxs;
	std::vector<std::array<Real,2>>                      midE_grid;
	std::vector<Real>                                    dcenter;
	std::vector<std::array<Integer,5>>                   C4E;
	std::vector<std::array<Integer,nside>>	             E4C;
	std::vector<std::array<Real,10>>	                 wtUt;
	std::vector<std::array<Integer,10>>                  Eindex;
	std::vector<std::array<Real,2>>                      areaV;
	std::vector<Real>   	                             areaE;
	std::vector<Real>   	                             areaI;
	std::vector<std::array<Real,nside>>                  LevdUn;
	std::vector<std::array<Real,nside>>                  NoutdUn;
	std::vector<std::array<Integer,nside>>               EindexVort;
	std::vector<Real>                                    arctheta;
	std::vector<std::array<Integer,nside>>               Edc;
	std::vector<std::array<Real,nside>>                  RAiv4T;
	std::vector<std::array<Real,nside>>                  RAiv4H;
	std::vector<std::array<Integer,nside>>               iVtopkey;
	std::vector<std::array<std::array<Real,2>,2>>        V_sph_grid;
	std::vector<Real>     lat;
	std::vector<Real>     lon;
	std::vector<Real>     latE;
	std::vector<Real>     lonE;
	std::vector<Real>     lon_v;
	std::vector<Real>     lat_v;
	std::vector<Real>     hb_a;
	std::vector<Real>     fcori;
	*/
	/*
	std::vector<Real>     icos_grid;
	std::vector<Real>     wtsph;
	std::vector<Real>     slen;
	std::vector<Real>     area;
	std::vector<Real>     Nvec;
	std::vector<Real>     Tvec;
	std::vector<Real>     Rvec;
	std::vector<Real>     Rcvec;
	std::vector<Integer>  nprox;
	std::vector<Integer>  prox;
	std::vector<Integer>  proxs;
	std::vector<Real>     midE_grid;
	std::vector<Real>     dcenter;
	std::vector<Integer>  C4E;
	std::vector<Integer>  E4C;
	std::vector<Real>     wtUt;
	std::vector<Integer>  Eindex;
	std::vector<Real>     areaV;
	std::vector<Real>     areaE;
	std::vector<Real>     areaI;
	std::vector<Real>     LevdUn;
	std::vector<Real>     NoutdUn;
	std::vector<Integer>  EindexVort;
	std::vector<Real>     arctheta;
	std::vector<Integer>  Edc;
	std::vector<Real>     RAiv4T;
	std::vector<Real>     RAiv4H;
	std::vector<Integer>  iVtopkey;
	std::vector<Real>     V_sph_grid;
	std::vector<Real>     lat;
	std::vector<Real>     lon;
	std::vector<Real>     latE;
	std::vector<Real>     lonE;
	std::vector<Real>     lon_v;
	std::vector<Real>     lat_v;
	std::vector<Real>     hb_a;
	std::vector<Real>     fcori;
	*/

	using int_array_1  = tbox::multi_array<Integer,1>;
	using int_array_2  = tbox::multi_array<Integer,2>;
	using int_array_3  = tbox::multi_array<Integer,3>;
	using real_array_1 = tbox::multi_array<Real,1>;
	using real_array_2 = tbox::multi_array<Real,2>;
	using real_array_3 = tbox::multi_array<Real,3>;

	static constexpr auto int_f_order_1  = int_array_1::f_order;
	static constexpr auto int_f_order_2  = int_array_2::f_order;
	static constexpr auto int_f_order_3  = int_array_3::f_order;
	static constexpr auto real_f_order_1 = real_array_1::f_order;
	static constexpr auto real_f_order_2 = real_array_2::f_order;
	static constexpr auto real_f_order_3 = real_array_3::f_order;

	real_array_2	icos_grid;
	real_array_3	wtsph;
	real_array_2	slen;
	real_array_1    area;
	real_array_3 	Nvec;
	real_array_3 	Tvec;
	real_array_3 	Rvec;
	real_array_2    Rcvec;
	int_array_1     nprox;
	int_array_2     prox;
	int_array_2     proxs;
	real_array_2    midE_grid;
	real_array_1    dcenter;
	int_array_2     C4E;
	int_array_2	    E4C;
	real_array_2	wtUt;
	int_array_2     Eindex;
	real_array_2    areaV;
	real_array_1   	areaE;
	real_array_1   	areaI;
	real_array_2    LevdUn;
	real_array_2    NoutdUn;
	int_array_2     EindexVort;
	real_array_1    arctheta;
	int_array_2     Edc;
	real_array_2    RAiv4T;
	real_array_2    RAiv4H;
	int_array_2     iVtopkey;
	real_array_3    V_sph_grid;
	real_array_1    lat;
	real_array_1    lon;
	real_array_1    latE;
	real_array_1    lonE;
	real_array_1    lon_v;
	real_array_1    lat_v;
	real_array_1    hb_a;
	real_array_1    fcori;

	Integer glvl;  // # grid level
	Integer nip;   // # of icosahedral points
	Integer niE;   // number of Edges
	Integer ndimFV;
	Integer NPTS;
	Integer ims, ime;
	Integer ips, ipe;
	Integer imsx, imex;
	Integer ims_E, ime_E;
	char    stagger;
};




#endif /* SRC_IFV_GRID_HPP_ */
