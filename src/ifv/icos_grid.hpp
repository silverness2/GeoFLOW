/*
 * grid.hpp
 *
 *  Created on: Jun 6, 2019
 *      Author: bflynt
 */

#ifndef SRC_IFV_GRID_HPP_
#define SRC_IFV_GRID_HPP_


#include <cstdint>
#include <string>
#include <vector>


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

	char stagger;
	Integer glvl;  // # grid level
	Integer nip;   // # of icosahedral points
	Integer niE;   // number of Edges
	static constexpr Integer ndim  = 3;
	static constexpr Integer nside = 6;
	Integer ndimFV;
	Integer NPTS;
	Integer ims, ime;
	Integer ips, ipe;
	Integer imsx, imex;
	Integer ims_E, ime_E;
};




#endif /* SRC_IFV_GRID_HPP_ */
