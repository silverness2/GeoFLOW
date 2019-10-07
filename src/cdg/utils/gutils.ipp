//==================================================================================
// Module       : gutils.ipp
// Date         : 1/31/19 (DLR)
// Description  : GeoFLOW utilities namespace
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

namespace geoflow
{

//**********************************************************************************
//**********************************************************************************
// METHOD : dopdf1d
// DESC   :
//          Computes probability distribution function (pdf) of 
//          scalar quantity u, and outputs it to the file specified.
//
// ARGS   :
//          grid  : grid object
//          u     : scalar field to find pdf for
//          tmp   : tmp vector of length at least X, each
//                  of same length as u
//          fname : file name for pdf
//          nbins : number of bins
//          bfixdr: if TRUE, will use fmin/fmax as bounds for dynamic
//                  range; else it will compute them dynamically
//          fmin/
//          fmax  : lowest and hight dynamic range value. This is 
//                  read or used if bfixdr = TRUE. If bfixdr = FALSE, 
//                  then full dynamic range of u is found, and fmin/fmax
//                  are returned with these values
//          dolog : if TRUE, take log of quantity magnitude when creating bins.
//                  default is FALSE.
// RETURNS: none.
//**********************************************************************************
template<typename T>
void dopdf1d(GGrid &grid, const GTVector<T> &u, GTVector<GTVector<T>*> &tmp, const GString &fname, GINT nbins, GBOOL bfixdr, T fmin, T fmax, GBOOL dolog);
{



} // end of method dopdf1d




} // end, namespace

