/*
 * grid.cpp
 *
 *  Created on: Jun 6, 2019
 *      Author: bflynt
 */

#include "icos_grid.hpp"
#include "xstd/fstream.hpp"

#include <cmath>

IcosGrid::IcosGrid()
	: stagger('A'){
}

void
IcosGrid::load(const char stag, const Integer level, const std::string filename){

	// Allocate Arrays
	stagger = stag;
	glvl    = level;
	nip     = 10 * std::pow(std::pow(2,glvl),2) + 2; // # of icosahedral points
	niE     = 30 * std::pow(4,glvl);                 // number of Edges
	ims     = 0;
	ime     = nip;
	ims_E   = 0;
	ime_E   = niE;

	if( stagger == 'A' ){
		ndimFV = 4;         //  Four vector  h/Vx/Vy/Vz
		ips    = 0;
		ipe    = nip;
		NPTS   = nip;
		imsx   = ims;
		imex   = ime;
	}
	else if( stagger == 'C' ){
		ndimFV = 3;         //  Four vector  h/Un/Ut, Ut is a linear combination of Un
		ips    = 0;
		ipe    = niE;
		NPTS   = niE;
		imsx   = ims_E;
		imex   = ime_E;
	}

	icos_grid.resize(2*nip);
	wtsph.resize(ndim*nside*nip);
	slen.resize(nside*nip);
	area.resize(nip);
	Nvec.resize(ndim*nside*nip);
	Tvec.resize(ndim*nside*nip);
	Rvec.resize(ndim*nside*nip);
	Rcvec.resize(ndim*nip);
	nprox.resize(nip);
	prox.resize(nside*nip);
	proxs.resize(nside*nip);
	midE_grid.resize(2*niE);
	dcenter.resize(niE);
	C4E.resize(5*niE);
	E4C.resize(nside*nip);
	wtUt.resize(10*niE);
	Eindex.resize(10*niE);
	areaV.resize(2*niE);
	areaE.resize(niE);
	areaI.resize(nip);
	LevdUn.resize(nside*niE);
	NoutdUn.resize(nside*nip);
	EindexVort.resize(nside*niE);
	arctheta.resize(niE);
	Edc.resize(nside*nip);
	RAiv4T.resize(nside*niE);
	RAiv4H.resize(nside*nip);
	iVtopkey.resize(nside*nip);
	V_sph_grid.resize(2*2*niE);

	lat.resize(nip);
	lon.resize(nip);
	latE.resize(niE);
	lonE.resize(niE);
	hb_a.resize(nip,0);
	//lon_v.resize(nip or niE);
	//lat_v.resize(nip or niE);

	// Open File
	xstd::iofort stream;
	stream.open(filename, xstd::iofort::ios::read);
	if( not stream.is_open() ){
		assert(false);
	}
	stream.read_record(
			std::make_pair(icos_grid.size(), icos_grid.data()),
			std::make_pair(wtsph.size(), wtsph.data()),
			std::make_pair(slen.size(), slen.data()),
			std::make_pair(area.size(), area.data()),
			std::make_pair(Nvec.size(), Nvec.data()),
			std::make_pair(Tvec.size(), Tvec.data()),
			std::make_pair(Rvec.size(), Rvec.data()),
			std::make_pair(Rcvec.size(), Rcvec.data()),
			std::make_pair(nprox.size(), nprox.data()),
			std::make_pair(prox.size(), prox.data()),
			std::make_pair(proxs.size(), proxs.data()),
			std::make_pair(midE_grid.size(), midE_grid.data()),
			std::make_pair(dcenter.size(), dcenter.data()),
			std::make_pair(C4E.size(), C4E.data()),
			std::make_pair(E4C.size(), E4C.data()),
			std::make_pair(wtUt.size(), wtUt.data()),
			std::make_pair(Eindex.size(), Eindex.data()),
			std::make_pair(areaV.size(), areaV.data()),
			std::make_pair(areaE.size(), areaE.data()),
			std::make_pair(areaI.size(), areaI.data()),
			std::make_pair(LevdUn.size(), LevdUn.data()),
			std::make_pair(NoutdUn.size(), NoutdUn.data()),
			std::make_pair(EindexVort.size(), EindexVort.data()),
			std::make_pair(arctheta.size(), arctheta.data()),
			std::make_pair(Edc.size(), Edc.data()),
			std::make_pair(RAiv4T.size(), RAiv4T.data()),
			std::make_pair(RAiv4H.size(), RAiv4H.data()),
			std::make_pair(iVtopkey.size(), iVtopkey.data()),
			std::make_pair(V_sph_grid.size(), V_sph_grid.data())
	);
	stream.close();

	for(Integer i = 0; i < nip; ++i){
		lon[i] = icos_grid[2*i];
		lat[i] = icos_grid[2*i+1];
	}

	for(Integer i = 0; i < niE; ++i){
		lonE[i] = midE_grid[2*i];
		latE[i] = midE_grid[2*i+1];
	}

	if(stagger == 'A'){
		lat_v = lat;
		lon_v = lon;
	}
	else if(stagger == 'C'){
		lat_v = latE;
		lon_v = lonE;
	}

}
