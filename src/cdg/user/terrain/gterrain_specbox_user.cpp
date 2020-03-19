// Module       : gterrain_specbox.hpp
// Date         : 7/11/19 (DLR)
// Description  : Terrain specification function implementations provided
//                by user for box grids (2d and 3d).
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#include "gterrain_specbox_user.hpp"


namespace gterrain_specbox {




//**********************************************************************************
//**********************************************************************************
// METHOD : impl_gauss_range
// DESC   : Specify terrange as a 'Gaussian mountain range' with
//          info by ttree. Intended for box grids
// ARGS   : ptree: main prop tree 
//          grid : grid
//          sblk : data block name
//          utmp : tmp arrays
//          xb   : terrain vectors
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_gauss_range(const PropertyTree &ptree, GString sblk, GGrid &grid, State &utmp, State &xb)
{

  PropertyTree ttree     = ptree.getPropertyTree(sblk);
  PropertyTree boxptree  = ptree.getPropertyTree("grid_box");

  GSIZET nxy;
  GFTYPE dx, dy, eps;
  GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();
  GTPoint<GFTYPE>  P0(3);
  
  nxy = (*xnodes)[0].size();

  std::vector<GFTYPE> x0   = ttree.getArray<GFTYPE>("x0");
  if ( GDIM == 3 ) 
  std::vector<GFTYPE> y0   = ttree.getArray<GFTYPE>("y0");
  std::vector<GFTYPE> xsig = ttree.getArray<GFTYPE>("xsigma");
  if ( GDIM == 3 ) 
  std::vector<GFTYPE> ysig = ttree.getArray<GFTYPE>("ysigma");
  std::vector<GFTYPE> h0   = ttree.getArray<GFTYPE>("h0");

  std::vector<GFTYPE> xyz0 = boxptree.getArray<GFTYPE>("xyz0");
//std::vector<GFTYPE> dxyz = boxptree.getArray<GFTYPE>("delxyz");
  P0 = xyz0; 

  eps =  100*std::numeric_limits<GFTYPE>::epsilon();

  // Build terrain height vector:
#if defined(_G_IS2D)
cout << "impl_gauss_range: made it!" << endl;
  for ( auto m=0; m<x0.size(); m++ ) {   // for each Gaussian lump
    *xb[1] = 0.0;
    for ( auto j=0; j<nxy; j++ ) {
      dx        = (*xnodes)[0][j] - x0[m];
      if ( FUZZYEQ(P0.x2,(*xnodes)[1][j],eps) ) {
        (*xb[1])[j] += h0[m]*exp(-pow(dx,2) / (xsig[m]*xsig[m]) );
cout << "impl_gauss_range: yb =" << (*xb[1])[j] << endl;
      }
    }
  }
#elif defined(_G_IS3D)
  *xb[2] = 0.0;
  for ( auto m=0; m<x0.size(); m++ ) {   // for each Gaussian lump
    for ( auto j=0; j<nxy; j++ ) {
      dx        = (*xnodes)[0][j] - x0[m];
      dy        = (*xnodes)[1][j] - y0[m];
      if ( FUZZYEQ(P0.x3,(*xnodes)[2][j],eps) ) {
        (*xb[2])[j] += h0[m]*exp(-dx*dx / (xsig[m]*xsig[m])) * exp(-dy*dy / (ysig[m]*ysig[m]));
      }
    }
  }
#endif
cout << "impl_gauss_range: xb =" << *xb[0] << endl;
cout << "impl_gauss_range: yb =" << *xb[1] << endl;

  return TRUE;

} // end, impl_gauss_range


//**********************************************************************************
//**********************************************************************************
// METHOD : impl_poly_range
// DESC   : Specify terrange as a inverse polynomials of the form:
//               h = h0/[(x-x0)^2/a^2 + (y-y0)^2/b^2]^(p/2)
//          
// ARGS   : ptree: main prop tree 
//          sblk : data block name
//          grid : grid
//          utmp : tmp arrays
//          xb   : terrain vectors
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_poly_range(const PropertyTree &ptree, GString sblk, GGrid &grid, State &utmp, State &xb)
{

  PropertyTree ttree     = ptree.getPropertyTree(sblk);
  PropertyTree boxptree  = ptree.getPropertyTree("grid_box");

  GSIZET nxy;
  GFTYPE dx, dy, eps;
  GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();
  GTPoint<GFTYPE>  P0(3);
  
  nxy = (*xnodes)[0].size();

  std::vector<GFTYPE> pexp = ttree.getArray<GFTYPE>("exponent");
  std::vector<GFTYPE> x0   = ttree.getArray<GFTYPE>("x0");
  if ( GDIM == 3 ) 
  std::vector<GFTYPE> y0   = ttree.getArray<GFTYPE>("y0");
  std::vector<GFTYPE> xsig = ttree.getArray<GFTYPE>("xsigma");
  if ( GDIM == 3 ) 
  std::vector<GFTYPE> ysig = ttree.getArray<GFTYPE>("ysigma");
  std::vector<GFTYPE> h0   = ttree.getArray<GFTYPE>("h0");

  std::vector<GFTYPE> xyz0 = boxptree.getArray<GFTYPE>("xyz0");
//std::vector<GFTYPE> dxyz = boxptree.getArray<GFTYPE>("delxyz");
  P0 = xyz0; 

  eps =  100*std::numeric_limits<GFTYPE>::epsilon();


  // Build terrain height vector:
#if defined(_G_IS2D)
  for ( auto m=0; m<x0.size(); m++ ) {   // for each Gaussian lump
    *xb[1] = 0.0;
    for ( auto j=0; j<nxy; j++ ) {
      dx        = (*xnodes)[0][j] - x0[m];
      if ( FUZZYEQ(P0.x2,(*xnodes)[1][j],eps) )
        (*xb[1])[j] += h0[m]/( pow( dx*dx/(xsig[m]*xsig[m]) + 1, pexp[m]/2) );
    }
  }
#elif defined(_G_IS3D)
  *xb[2] = 0.0;
  for ( auto m=0; m<x0.size(); m++ ) {   // for each Gaussian lump
    for ( auto j=0; j<nxy; j++ ) {
      dx        = (*xnodes)[0][j] - x0[m];
      dy        = (*xnodes)[1][j] - y0[m];
      if ( FUZZYEQ(P0.x3,(*xnodes)[2][j],eps) )
        (*xb[2])[j] += h0[m]/( pow( dx*dx/(xsig[m]*xsig[m]) + dy*dy/(ysig[m]*ysig[m]) + 1, pexp[m]/2) );
    }
  }
#endif

  return TRUE;

} // end, impl_poly_range



//**********************************************************************************
//**********************************************************************************
// METHOD : impl_schar_range
// DESC   : Specify terrange as given in Schar et al., MWR 130:2459 2002:
//               h(x) = cos^2(pi x/lambda) h'(x),
//          where
//               h'(x) = h0 cos^2(pi x/2a), |x| <= a;
//                     = 0                , |x| >  a.
//          Currently works only for 2D.
//          
// ARGS   : ptree: main prop tree 
//          sblk : data block name
//          grid : grid
//          utmp : tmp arrays
//          xb   : terrain vectors
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_schar_range(const PropertyTree &ptree, GString sblk, GGrid &grid, State &utmp, State &xb)
{


  PropertyTree ttree     = ptree.getPropertyTree(sblk);
  PropertyTree boxptree  = ptree.getPropertyTree("grid_box");

  GSIZET nxy;
  GFTYPE eps, x;
  GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();
  GTPoint<GFTYPE>  P0(3);
  
  nxy = (*xnodes)[0].size();

  GFTYPE lambda  = ttree.getValue<GFTYPE>("lambda");       // perturbation wavelength
  GFTYPE extent  = ttree.getValue<GFTYPE>("range_extent"); // extent 
  GFTYPE h0      = ttree.getValue<GFTYPE>("h0");           // height
  std::vector<GFTYPE> xyz0 = boxptree.getArray<GFTYPE>("xyz0");
//std::vector<GFTYPE> dxyz = boxptree.getArray<GFTYPE>("delxyz");
  P0 = xyz0; 

  eps =  100*std::numeric_limits<GFTYPE>::epsilon();


  // Build terrain height vector:
#if defined(_G_IS2D)
  *xb[1] = 0.0;
  for ( auto j=0; j<nxy; j++ ) {
    x        = (*xnodes)[0][j];
    if ( FUZZYEQ(P0.x2,(*xnodes)[1][j],eps) )
//               h(x) = cos^2(pi x/lambda) h'(x),
//               h'(x) = h0 cos^2(pi x/2a), |x| <= a;
      (*xb[1])[j] = 0.0;
      if ( abs(x) <= extent ) {
        (*xb[1])[j] += h0*pow(cos(PI*x/lambda)    ,2)
                         *pow(cos(PI*x/(2*extent)),2);
      }
  }
#elif defined(_G_IS3D)
  assert(FALSE && "Method undefined in 3D");
#endif

  return TRUE;

} // end, impl_schar_range




} // end, ginitbdy namespace
