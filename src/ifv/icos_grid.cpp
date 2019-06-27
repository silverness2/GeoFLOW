/*
 * grid.cpp
 *
 *  Created on: Jun 6, 2019
 *      Author: bflynt
 */

//#include "tbox/pio.hpp"


#include "icos_grid.hpp"
#include "xstd/fstream.hpp"


#include <cmath>

IcosGrid::IcosGrid() {
}

void
IcosGrid::load(const char stag, const Integer level, const std::string filename){

	// Allocate Arrays
	stagger = stag;
	glvl    = level;
	nip     = 10 * std::pow(std::pow(2,glvl),2) + 2; // # of icosahedral points
	niE     = 30 * std::pow(4,glvl);                 // number of Edges
	ims     = 1;
	ime     = nip;
	ims_E   = 1;
	ime_E   = niE;

	if( stagger == 'A' ){
		ndimFV = 4;         //  Four vector  h/Vx/Vy/Vz
		ips    = 1;
		ipe    = nip;
		NPTS   = nip;
		imsx   = ims;
		imex   = ime;
	}
	else if( stagger == 'C' ){
		ndimFV = 3;         //  Four vector  h/Un/Ut, Ut is a linear combination of Un
		ips    = 1;
		ipe    = niE;
		NPTS   = niE;
		imsx   = ims_E;
		imex   = ime_E;
	}

	/*
	icos_grid.resize(nip);        // icos_grid.resize(2*nip);
	wtsph.resize(nip);            // wtsph.resize(ndim*nside*nip);
	slen.resize(nip);             // slen.resize(nside*nip);
	area.resize(nip);
	Nvec.resize(nip);             // Nvec.resize(ndim*nside*nip);
	Tvec.resize(nip);             // Tvec.resize(ndim*nside*nip);
	Rvec.resize(nip);             // Rvec.resize(ndim*nside*nip);
	Rcvec.resize(nip);            // Rcvec.resize(ndim*nip);
	nprox.resize(nip);
	prox.resize(nip);             // prox.resize(nside*nip);
	proxs.resize(nip);            // proxs.resize(nside*nip);
	midE_grid.resize(niE);        // midE_grid.resize(2*niE);
	dcenter.resize(niE);
	C4E.resize(niE);              // C4E.resize(5*niE);
	E4C.resize(nip);              // E4C.resize(nside*nip);
	wtUt.resize(niE);             // wtUt.resize(10*niE);
	Eindex.resize(niE);           // Eindex.resize(10*niE);
	areaV.resize(niE);            // areaV.resize(2*niE);
	areaE.resize(niE);
	areaI.resize(nip);
	LevdUn.resize(niE);           // LevdUn.resize(nside*niE);
	NoutdUn.resize(nip);          // NoutdUn.resize(nside*nip);
	EindexVort.resize(niE);       // EindexVort.resize(nside*niE);
	arctheta.resize(niE);
	Edc.resize(nip);              // Edc.resize(nside*nip);
	RAiv4T.resize(niE);           // RAiv4T.resize(nside*niE);
	RAiv4H.resize(nip);           // RAiv4H.resize(nside*nip);
	iVtopkey.resize(nip);         // iVtopkey.resize(nside*nip);
	V_sph_grid.resize(niE);       // V_sph_grid.resize(2*2*niE);
	 */
	/*
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
	*/

	icos_grid.resize(2,nip).reorder(real_f_order_2).rebase(1);
	wtsph.resize(ndim,nside,nip).reorder(real_f_order_3).rebase(1);
	slen.resize(nside,nip).reorder(real_f_order_2).rebase(1);
	area.resize(nip).reorder(real_f_order_1).rebase(1);
	Nvec.resize(ndim,nside,nip).reorder(real_f_order_3).rebase(1);
	Tvec.resize(ndim,nside,nip).reorder(real_f_order_3).rebase(1);
	Rvec.resize(ndim,nside,nip).reorder(real_f_order_3).rebase(1);
	Rcvec.resize(ndim,nip).reorder(real_f_order_2).rebase(1);
	nprox.resize(nip).reorder(int_f_order_1).rebase(1);
	prox.resize(nside,nip).reorder(int_f_order_2).rebase(1);
	proxs.resize(nside,nip).reorder(int_f_order_2).rebase(1);
	midE_grid.resize(2,niE).reorder(real_f_order_2).rebase(1);
	dcenter.resize(niE).reorder(real_f_order_1).rebase(1);
	C4E.resize(5,niE).reorder(int_f_order_2).rebase(1);
	E4C.resize(nside,nip).reorder(int_f_order_2).rebase(1);
	wtUt.resize(10,niE).reorder(real_f_order_2).rebase(1);
	Eindex.resize(10,niE).reorder(int_f_order_2).rebase(1);
	areaV.resize(2,niE).reorder(real_f_order_2).rebase(1);
	areaE.resize(niE).reorder(real_f_order_1).rebase(1);
	areaI.resize(nip).reorder(real_f_order_1).rebase(1);
	LevdUn.resize(nside,niE).reorder(real_f_order_2).rebase(1);
	NoutdUn.resize(nside,nip).reorder(real_f_order_2).rebase(1);
	EindexVort.resize(nside,niE).reorder(int_f_order_2).rebase(1);
	arctheta.resize(niE).reorder(real_f_order_1).rebase(1);
	Edc.resize(nside,nip).reorder(int_f_order_2).rebase(1);
	RAiv4T.resize(nside,niE).reorder(real_f_order_2).rebase(1);
	RAiv4H.resize(nside,nip).reorder(real_f_order_2).rebase(1);
	iVtopkey.resize(nside,nip).reorder(int_f_order_2).rebase(1);
	V_sph_grid.resize(2,2,niE).reorder(real_f_order_3).rebase(1);

	lat.resize(nip).reorder(real_f_order_1).rebase(1);
	lon.resize(nip).reorder(real_f_order_1).rebase(1);
	latE.resize(niE).reorder(real_f_order_1).rebase(1);
	lonE.resize(niE).reorder(real_f_order_1).rebase(1);
	hb_a.resize(nip).reorder(real_f_order_1).rebase(1).fill(0);
	fcori.resize(NPTS).reorder(real_f_order_1).rebase(1).fill(0);
	//lon_v.resize(nip or niE);
	//lat_v.resize(nip or niE);

	// Open File
	xstd::iofort stream;
	stream.open(filename, xstd::iofort::ios::read);
	if( not stream.is_open() ){
		assert(false);
	}
	stream.read_record(
			std::make_pair( icos_grid.size(), icos_grid.data()),
			std::make_pair( wtsph.size(), wtsph.data()),
			std::make_pair( slen.size(), slen.data()),
			std::make_pair( area.size(), area.data()),
			std::make_pair( Nvec.size(), Nvec.data()),
			std::make_pair( Tvec.size(), Tvec.data()),
			std::make_pair( Rvec.size(), Rvec.data()),
			std::make_pair( Rcvec.size(), Rcvec.data()),
			std::make_pair( nprox.size(), nprox.data()),
			std::make_pair( prox.size(), prox.data()),
			std::make_pair( proxs.size(), proxs.data()),
			std::make_pair( midE_grid.size(), midE_grid.data()),
			std::make_pair( dcenter.size(), dcenter.data()),
			std::make_pair( C4E.size(), C4E.data()),
			std::make_pair( E4C.size(), E4C.data()),
			std::make_pair( wtUt.size(), wtUt.data()),
			std::make_pair( Eindex.size(), Eindex.data()),
			std::make_pair( areaV.size(), areaV.data()),
			std::make_pair( areaE.size(), areaE.data()),
			std::make_pair( areaI.size(), areaI.data()),
			std::make_pair( LevdUn.size(), LevdUn.data()),
			std::make_pair( NoutdUn.size(), NoutdUn.data()),
			std::make_pair( EindexVort.size(), EindexVort.data()),
			std::make_pair( arctheta.size(), arctheta.data()),
			std::make_pair( Edc.size(), Edc.data()),
			std::make_pair( RAiv4T.size(), RAiv4T.data()),
			std::make_pair( RAiv4H.size(), RAiv4H.data()),
			std::make_pair( iVtopkey.size(), iVtopkey.data()),
			std::make_pair( V_sph_grid.size(), V_sph_grid.data())
	);
	stream.close();

	//using namespace geoflow::tbox;
	//pio::pout << "icos_grid.data() = " << std::endl;
	//pio::pout << icos_grid.data()[0] << std::endl;
	//pio::pout << icos_grid.data()[1] << std::endl;
	//pio::pout << icos_grid.data()[2] << std::endl;
	//pio::pout << icos_grid.data()[3] << std::endl;


	for(Integer ip = 1; ip <= nip; ++ip){
		lon(ip) = icos_grid(1,ip);
		lat(ip) = icos_grid(2,ip);
	}

	for(Integer iE = 1; iE <= niE; ++iE){
        lonE(iE) = midE_grid(1,iE);
        latE(iE) = midE_grid(2,iE);
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
