//==================================================================================
// Module       : ginitstate_user.cpp
// Date         : 7/10/19 (DLR)
// Description  : Direct user state initialization function implementations. These
//                methods are called directly during configuration, and can set 
//                forcing for entire state (v+b+s) etc. The 'component' types
//                in which component groups (v, b, s, etc) are set individually
//                are contained in separate namespaces.
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#include "gmtk.hpp"
#include "ginitstate_direct_user.hpp"


namespace ginitstate {



//**********************************************************************************
//**********************************************************************************
// METHOD : impl_boxnwaveburgers
// DESC   : Initialize state for Burgers with N-wave on box grids with
//          Dirichlet or periodic boundaries
// ARGS   : ptree  : main prop tree
//          sconfig: ptree block name containing variable config
//          grid   : grid
//          t      : time
//          utmp   : tmp arrays
//          ub     : bdy vectors (one for each state element)
//          u      : current state
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_boxnwaveburgers(const PropertyTree &ptree, GString &sconfig, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{
  GString          serr = "impl_boxnwaveburgers: ";
  GBOOL            bret = TRUE, brot = FALSE;
  GINT             nlump=0;
  GSIZET           i, j, nxy;
  GFTYPE           K2, nu, Re, r2, tdenom;
  GFTYPE           efact, sum, tfact, xfact;
  GTVector<GFTYPE> xx(GDIM), si(GDIM), sig(GDIM), t0;
  GTPoint<GFTYPE>  kprop(3), r0(3), P0(3), gL(3);
  std::vector<GFTYPE>  kxprop, kyprop, kzprop;
  std::vector<GFTYPE>  xinit , yinit , zinit ;
  std::vector<GBOOL>   bplanar;
  std::vector<GFTYPE>  tinit;
  std::vector<GFTYPE>  ULparam;
  GString              snut;

  PropertyTree nwaveptree = ptree   .getPropertyTree(sconfig);
  PropertyTree boxptree   = ptree   .getPropertyTree("grid_box");
  PropertyTree nuptree    = ptree.getPropertyTree("dissipation_traits");

  GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();

  assert(grid.gtype() == GE_REGULAR && "Invalid element types");

  GTVector<GString> bc(6);
  bc[0] = boxptree.getValue<GString>("bdy_x_0");
  bc[1] = boxptree.getValue<GString>("bdy_x_1");
  bc[2] = boxptree.getValue<GString>("bdy_y_0");
  bc[3] = boxptree.getValue<GString>("bdy_y_1");
  bc[4] = boxptree.getValue<GString>("bdy_z_0");
//    && "INFLOWT boundaries must be set on all boundaries");

  nxy = (*xnodes)[0].size(); // same size for x, y, z

  // From Whitham's book, in 1d:
  // u(x,t) = (x/t) [ 1 + sqrt(t/t0) (e^Re - 1)^-1 exp(x^2/(4 nu t))i ]^-1
  // were Re is 'Reynolds' number: Re = A / 2nu; can think of
  // A ~ U L scaling. But we won't parameterize in terms of Re, 
  // but rather, nu.
  // Get some parameters; xinit, tinit, ULparam, bplanar, kprop,
  // should have the same number, nlump, elements, one foreach 'wave':
  xinit      = nwaveptree.getArray<GFTYPE>("x0"); 
  yinit      = nwaveptree.getArray<GFTYPE>("y0"); 
  zinit      = nwaveptree.getArray<GFTYPE>("z0"); 
  tinit      = nwaveptree.getArray<GFTYPE>("t0");
  ULparam    = nwaveptree.getArray<GFTYPE>("ULparam");
  bplanar    = nwaveptree.getArray<GBOOL> ("planar");
  kxprop = nwaveptree.getArray<GFTYPE>("prop_dir_x");
  kyprop = nwaveptree.getArray<GFTYPE>("prop_dir_y");
  kzprop = nwaveptree.getArray<GFTYPE>("prop_dir_z");
//Re         = nwaveptree.getValue<GFTYPE>("Re",6.0);

  nu       = nuptree   .getValue<GFTYPE>("nu",0.0833);
  snut     = nuptree   .getValue<GString>("nu_type","constant");

  t0 = tinit;

  assert(yinit  .size() == xinit.size()
      && tinit  .size() == xinit.size()
      && ULparam.size() == xinit.size()
      && bplanar.size() == xinit.size()
      && kxprop .size() == xinit.size()
      && kyprop .size() == xinit.size()
      && "(1)Lump count must be consistent");

  if ( GDIM > 2 ) {
    assert(zinit  .size() == xinit.size()
        && kzprop .size() == xinit.size()
        && "(2)Lump count must be consistent");
  }
  nlump = xinit.size();

  assert( snut == "constant" && "nu_type must bet set to 'constant')");
  cout << "impl_boxnwaveburgers: nu=" << nu << " Re=" << Re << " time=" << time << endl;

  if ( time <= 10.0*std::numeric_limits<GFTYPE>::epsilon() ) time = t0.max();

  for ( i=0; i<GDIM; i++ ) {
    *u[i] = 0.0;
  }

  for ( GINT ilump=0; ilump<nlump; ilump++ ) {
    r0[0]  = xinit[ilump]; r0[1]  = yinit[ilump]; 
    if ( GDIM > 2 ) r0[2]  = zinit[ilump]; 
    kprop[0] = kxprop[ilump]; kprop[1] = kyprop[ilump];
    if ( GDIM > 2 ) kprop[2]  = kzprop[ilump]; 
    Re  = ULparam[ilump]/nu; // set Re from nu and ULparam
    kprop  *= 1.0/kprop.norm();
    tdenom  = 1.0/(4.0*nu*time);
    tfact   = bplanar[ilump] ? sqrt(time/t0[ilump]): time/t0[ilump];
    for ( i=0, K2=0.0; i<GDIM; i++ ) K2 += kprop[i]*kprop[i];

    // If prop direction has more than one component != 0. Then
    // front is rotated (but still planar):
//  for ( i=0, brot=TRUE; i<GDIM; i++ ) brot = brot && K[i] != 0.0 ;
//  for ( i=0, idir=0; i<GDIM; i++ ) if ( K[i] > 0 ) {idir=i; break;}
//  K2 = brot && K2 == 0.0 ? 1.0 : K2;
//  if ( time <= 10.0*std::numeric_limits<GFTYPE>::epsilon() ) time = K2 * t0;
    for ( j=0; j<nxy; j++ ) {
      for ( i=0; i<GDIM; i++ ) {
        xx[i] = (*xnodes)[i][j] - r0[i];
      }
      if ( bplanar[ilump] ) { // compute k.r for planar wave
        for ( i=0, sum=0.0; i<GDIM; i++ ) { 
          sum += kprop[i]*xx[i];
          xx[i] = 0.0;
        }
        xx[0] = sum;
      }
      for ( i=0, r2=0.0; i<GDIM; i++ ) r2 += xx[i]*xx[i];  
  
      efact   = tfact * exp(r2*tdenom) / ( exp(Re) - 1.0 );
      xfact   = 1.0 /( time * (  1.0 + efact ) );
      // u(x,t) = (x/t) [ 1 + sqrt(t/t0) (e^Re - 1)^-1 exp(x^2/(4 nu t)) ]^-1
      for ( i=0; i<GDIM; i++ ) (*u[i])[j] += xx[i]*xfact;
  //cout << "impl_boxnwaveburgers: ux[" << j << "]=" << (*u[0])[j] << endl;
      // dU1max = 1.0 / ( time * (sqrt(time/A) + 1.0) );
      // aArea  = 4.0*nu*log( 1.0 + sqrt(A/time) );
    } // end, coord loop
  } // end, lump loop

  return bret;

} // end, impl_icosnwaveburgers


//**********************************************************************************
//**********************************************************************************
// METHOD : impl_icosnwaveburgers
// DESC   : Initialize state for Burgers with N-wave on icos grids
// ARGS   : ptree  : main prop tree
//          sconfig: ptree block name containing variable config
//          grid   : grid
//          t      : time
//          utmp   : tmp arrays
//          ub     : bdy vectors (one for each state element)
//          u      : current state
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_icosnwaveburgers(const PropertyTree &ptree, GString &sconfig, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{
  GString          serr = "impl_icosnwaveburgers: ";
  GBOOL            bret;
  GSIZET           i, j, nxy;
  GFTYPE           nu, Re, r, s, tdenom;
  GFTYPE           lat, lon;
  GFTYPE           efact, sum, tfact, xfact;
  GFTYPE           x, y, z;
  GTVector<GFTYPE>            t0, xx(GDIM+1);
  GTVector<GTPoint<GFTYPE>>   r0(GDIM+1);
  std::vector<GFTYPE>         lat0, lon0, st0, Uparam;

  PropertyTree nwaveptree = ptree   .getPropertyTree(sconfig);
  PropertyTree gridptree  = ptree   .getPropertyTree("grid_icos");
  PropertyTree nuptree    = ptree.getPropertyTree("dissipation_traits");

  GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();

  assert(grid.gtype() == GE_2DEMBEDDED && "Invalid element types");
  assert(u.size() >= GDIM+1 && "Insufficient number of state members");


  nxy = (*xnodes)[0].size(); // same size for x, y, z

  // From Whitham's book, in 1d:
  // u(x,t) = (x/t) [ 1 + sqrt(t/t0) (e^Re - 1)^-1 exp(x^2/(4 nu t))i ]^-1
  // were Re is 'Reynolds' number: Re = A / 2nu; can think of
  // A ~ U L scaling. But we won't parameterize in terms of Re, 
  // but rather, nu.
  // Set some parameters:
  r      = gridptree.getValue <GFTYPE>("radius");
  lat0   = nwaveptree.getArray<GFTYPE>("latitude0"); 
  lon0   = nwaveptree.getArray<GFTYPE>("longitude0"); 
  Uparam = nwaveptree.getArray<GFTYPE>("Uparam");
//Re     = nwaveptree.getValue<GFTYPE>("Re",6.0);
  st0    = nwaveptree.getArray<GFTYPE>("t0");
  nu     = nuptree   .getValue<GFTYPE>("nu",0.0833);

  t0.resize(st0.size());
  t0     = st0;
  assert(lat0.size() == lon0.size() 
      && lat0.size() == t0.size()
      && "lat0, lon0, and t0 must be the same size");

  if ( time <= 10.0*std::numeric_limits<GFTYPE>::epsilon() ) time = t0.max();

  // Convert initial positions to radians:
  for ( GINT ilump=0; ilump<lat0.size(); ilump++) {
    lat0[ilump]    *= (PI/180.0);
    lon0[ilump]    *= (PI/180.0);
    r0  [ilump].x1  = r*cos(lat0[ilump])*cos(lon0[ilump]);
    r0  [ilump].x2  = r*cos(lat0[ilump])*sin(lon0[ilump]);
    r0  [ilump].x3  = r*sin(lat0[ilump]);
  }

  for ( i=0; i<GDIM+1; i++ ) *u[i] = 0.0;
   
  // Initialize each lump:
  for ( GINT ilump=0; ilump<lat0.size(); ilump++) {
    Re = Uparam[ilump]*r/nu; // set Re from nu, U, radius
    for ( j=0; j<nxy; j++ ) {
       x = (*xnodes)[0][j]; y = (*xnodes)[1][j]; z = (*xnodes)[2][j];
       lat = asin(z/r);
       lon = atan2(y,x);
      for ( i=0; i<GDIM+1; i++ ) { 
        xx[i] = (*xnodes)[i][j] - r0[ilump][i];
      }
/*
      xx[0] = r*lat;
      xx[1] = r*lon;
*/
      // find arclength from lump center
      s     = r*acos( sin(lat)*sin(lat0[ilump])
            +       cos(lat)*cos(lat0[ilump])*cos(lon - lon0[ilump]) );
  
      tdenom = 1.0/(4.0*nu*time);
      tfact  = time/t0[ilump];
      efact  = tfact * exp(s*s*tdenom) / ( exp(Re) - 1.0 );
      xfact  = 1.0 /( time * (  1.0 + efact ) );
      for ( i=0; i<GDIM+1; i++ ) {
        (*u[i])[j] += xx[i]*xfact;
        assert( std::isfinite((*u[i])[j]) );
      }
      // dU1max = 1.0 / ( time * (sqrt(time/A) + 1.0) );
      // aArea  = 4.0*nu*log( 1.0 + sqrt(A/time) );
    } // end, coord j-loop 
  } // end, ilump-loop

//GMTK::vsphere2cart(grid, usph, GVECTYPE_PHYS, u);
  GMTK::constrain2sphere(grid, u);

  bret = TRUE;
  for ( j=0; j<GDIM+1; j++ ) {
     bret = bret && u[j]->isfinite(i);
  }

  assert(bret && "Initial conditions not finite!");


  return bret;

} // end, impl_icosnwaveburgers



//**********************************************************************************
//**********************************************************************************
// METHOD : impl_boxdirgauss
// DESC   : Initialize state for Burgers with Gauss lump on box grids with
//          Dirichlet boundaries
// ARGS   : ptree  : main prop tree
//          sconfig: ptree block name containing variable config
//          grid   : grid
//          time   : time
//          utmp   : tmp arrays
//          ub     : bdy vectors (one for each state element)
//          u      : current state
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_boxdirgauss(const PropertyTree &ptree, GString &sconfig, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{
  GString          serr = "impl_boxdirgauss: ";
  GBOOL            bContin;
  GINT             j, n;
  GFTYPE           argxp;
  GFTYPE           nxy, nu, sig0, E0;
  GTVector<GFTYPE> xx(GDIM), si(GDIM), sig(GDIM), ufact(GDIM);
  State            c(GDIM);
  GTPoint<GFTYPE>  r0(3), P0(3);
  GString          snut;

  PropertyTree heatptree = ptree.getPropertyTree("sconfig");
  PropertyTree boxptree  = ptree.getPropertyTree("grid_box");
  PropertyTree advptree  = ptree.getPropertyTree("pde_burgers");
  PropertyTree nuptree   = ptree.getPropertyTree("dissipation_traits");
  GBOOL doheat   = advptree.getValue<GBOOL>("doheat");
  GBOOL bpureadv = advptree.getValue<GBOOL>("bpureadv");

  assert((bpureadv || doheat) && "Pure advection or heat must be used");

  GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();

  assert(grid.gtype() == GE_REGULAR && "Invalid element types");

  std::vector<GFTYPE> cs;
  if ( bpureadv ) {
    cs = heatptree.getArray<GFTYPE>("adv_vel");
  }

  for ( GSIZET j=0; j<GDIM; j++ ) c[j] = u[j+1];

  // Check bdy conditioins:
  GTVector<GString> bc(6);
  bc[0] = boxptree.getValue<GString>("bdy_x_0");
  bc[1] = boxptree.getValue<GString>("bdy_x_1");
  bc[2] = boxptree.getValue<GString>("bdy_y_0");
  bc[3] = boxptree.getValue<GString>("bdy_y_1");
  bc[4] = boxptree.getValue<GString>("bdy_z_0");
  bc[5] = boxptree.getValue<GString>("bdy_z_1");
  assert(bc.multiplicity("GBDY_INFLOWT") >= 2*GDIM
      && "GBDY_INFLOW boundaries must be set on all boundaries");

  nxy = (*xnodes)[0].size(); // same size for x, y, z

  r0.x1 = heatptree.getValue<GFTYPE>("x0");
  r0.x2 = heatptree.getValue<GFTYPE>("y0");
  r0.x3 = heatptree.getValue<GFTYPE>("z0");
  sig0  = heatptree.getValue<GFTYPE>("sigma");
  E0    = heatptree.getValue<GFTYPE>("E0");

  nu     = nuptree   .getValue<GFTYPE>("nu");
  snut   = nuptree   .getValue<GString>("nu_type","constant");
  assert( snut == "constant" && "nu_type must bet set to 'constant')");

  // Set velocity here. May be a function of time.
  // These point to components of state u_:
  for ( j=0; j<GDIM; j++ ) *c[j] = 0.0;

  if ( bpureadv ) {
     for ( j=0; j<GDIM; j++ ) *c[j] = cs[j];
  }


  // Prepare for case where sig is anisotropic (for later, maybe):
  for ( GSIZET k=0; k<GDIM; k++ ) {
    sig  [k] = sqrt(sig0*sig0 + 4.0*time*nu); // constant viscosity only
    si   [k] = 1.0/(sig[k]*sig[k]);
    ufact[k] = sqrt(2*E0)*pow(sig0/sig[k],GDIM);
  }

  // Ok, return to assumption of isotropic nu: 
  for ( GSIZET j=0; j<nxy; j++ ) {
    // Note: following c t is actually Integral_0^t c(t') dt', 
    //       so if c(t) changes, change this term accordingly:
    for ( GSIZET i=0; i<GDIM; i++ ) xx[i] = (*xnodes)[i][j] - r0[i] - (*c[i])[j]*time;
    argxp = 0.0;
    for ( GSIZET i=0; i<GDIM; i++ ) argxp += -pow(xx[i],2.0)*si[i];
   (*u[0])[j] = ufact[0]*exp(argxp);
  }

  return TRUE;
} // end, impl_boxdirgauss


//**********************************************************************************
//**********************************************************************************
// METHOD : impl_boxpergauss
// DESC   : Initialize state for Burgers with Gauss lump on box grids with
//          periodic boundaries
// ARGS   : stree  : main prop tree
//          sconfig: ptree block name containing variable config
//          grid   : grid
//          t      : time
//          utmp   : tmp arrays
//          ub     : bdy vectors (one for each state element)
//          u      : current state
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_boxpergauss(const PropertyTree &ptree, GString &sconfig, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{
  GString          serr = "impl_boxpergauss: ";
  GBOOL            bContin;
  GSIZET           i, j, k, n;
  GFTYPE           iargp, iargm ;
  GFTYPE           isum , irat , prod;
  GFTYPE           sumn , eps;
  GFTYPE           nxy, nu, pint, sig0, E0;
  GTVector<GFTYPE> f(GDIM), xx(GDIM), si(GDIM), sig(GDIM);
  GTPoint<GFTYPE>  r0(3), P0(3), gL(3);
  State            c(GDIM);
  GString          snut;

  PropertyTree heatptree = ptree.getPropertyTree(sconfig);
  PropertyTree boxptree  = ptree.getPropertyTree("grid_box");
  PropertyTree advptree  = ptree.getPropertyTree("pde_burgers");
  PropertyTree nuptree   = ptree.getPropertyTree("dissipation_traits");
  GBOOL doheat   = advptree.getValue<GBOOL>("doheat");
  GBOOL bpureadv = advptree.getValue<GBOOL>("bpureadv");
  
  assert((bpureadv || doheat) && "Pure advection or heat must be used");

  GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();

  assert(grid.gtype() == GE_REGULAR && "Invalid element types");

  eps = 1.0e-4*std::numeric_limits<GFTYPE>::epsilon();

  // Get periodicity length, gL:
  std::vector<GFTYPE> xyz0 = boxptree.getArray<GFTYPE>("xyz0");
  std::vector<GFTYPE> dxyz = boxptree.getArray<GFTYPE>("delxyz");
  P0 = xyz0; r0 = dxyz; gL = r0;

  for ( GSIZET j=0; j<GDIM; j++ ) c[j] = u[j+1];

  std::vector<GFTYPE> cs;
  if ( bpureadv ) {
    cs = heatptree.getArray<GFTYPE>("adv_vel");
  }

  GTVector<GString> bc(6);
  bc[0] = boxptree.getValue<GString>("bdy_x_0");
  bc[1] = boxptree.getValue<GString>("bdy_x_1");
  bc[2] = boxptree.getValue<GString>("bdy_y_0");
  bc[3] = boxptree.getValue<GString>("bdy_y_1");
  bc[4] = boxptree.getValue<GString>("bdy_z_0");
  bc[5] = boxptree.getValue<GString>("bdy_z_1");
  assert(bc.multiplicity("GBDY_PERIODIC") >= 2*GDIM
      && "Periodic boundaries must be set on all boundaries");

  nxy = (*xnodes)[0].size(); // same size for x, y, z

  r0.x1 = heatptree.getValue<GFTYPE>("x0");
  r0.x2 = heatptree.getValue<GFTYPE>("y0");
  r0.x3 = heatptree.getValue<GFTYPE>("z0");
  sig0  = heatptree.getValue<GFTYPE>("sigma");
  E0    = heatptree.getValue<GFTYPE>("E0");

  nu    = nuptree   .getValue<GFTYPE>("nu");
  snut  = nuptree   .getValue<GString>("nu_type","constant");
  assert( snut == "constant" && "nu_type must bet set to 'constant')");


  // Set adv velocity components. Note:
  // First state elem is the scalar solution, and
  // the remainder are the velocity components:

  for ( j=0; j<GDIM; j++ ) {
    sig[j] = sqrt(sig0*sig0 + 4.0*time*nu);
    si [j] = 1.0/(sig[j]*sig[j]);
   *c  [j] = 0.0;
  }

  // Set velocity here. May be a function of time.
  // These point to components of state u_:
  if ( bpureadv ) for ( j=0; j<GDIM; j++ ) *c[j] = cs[j];

  for ( n=0; n<nxy; n++ ) {

    prod = 1.0;
    for ( k=0; k<GDIM; k++ ) {
      // Note: following c t is actually Integral_0^t c(t') dt', 
      //       so if c(t) changes, change this term accordingly:
      f [k]  = modf((*c[k])[j]*time/gL[k],&pint);
//    f [k]  = (*c[k])[n]*t/gL[k];
      xx[k]  = (*xnodes)[k][n] - r0[k] - f[k]*gL[k];

      isum    = 0.0;
      i       = 0;
      irat    = 1.0;
      while ( irat > eps ) { // inner sum
        iargp   = pow((xx[k]+i*gL[k]),2)*si[k];
        iargm   = pow((xx[k]-i*gL[k]),2)*si[k];
        sumn    = i==0 ? exp(-iargp) : exp(-iargp) + exp(-iargm);
        isum   += sumn;
        irat    = sumn / isum ;
        i++;
      }
      prod *= isum;
    }
    (*u[0])[n] = sqrt(2.0*E0)*pow(sig0,GDIM)/pow(sig[0],GDIM)*prod;

  } // end, loop over grid points

  return TRUE;

} // end, impl_boxpergauss


//**********************************************************************************
//**********************************************************************************
// METHOD : impl_icosgauss
// DESC   : Initialize state for Burgers with Gauss lump on ICOS grid
// ARGS   : ptree  : main prop tree
//          sconfig: ptree block name containing variable config
//          grid   : grid
//          t      : time
//          utmp   : tmp arrays
//          ub     : bdy vectors (one for each state element)
//          u      : current state
// RETURNS: TRUE on success; else FALSE 
//**********************************************************************************
GBOOL impl_icosgauss(const PropertyTree &ptree, GString &sconfig, GGrid &grid, Time &time, State &utmp, State &ub, State &u)
{

  GString          serr = "impl_icosgauss: ";
  GBOOL            bContin;
  GINT             j, k, n, nlumps;
  GSIZET           nxy;
  GFTYPE           alpha, argxp;
  GFTYPE           lat, lon;
  GFTYPE           x, y, z, r, s;
  GFTYPE           nu, rad, E0, u0 ;
  GFTYPE           vtheta, vphi;
  GFTYPE           tiny = std::numeric_limits<GFTYPE>::epsilon();
  GTPoint<GFTYPE>           rt(3);
  GTVector<GFTYPE>          xx(3);
  GTVector<GFTYPE>          si(4), sig(4), ufact(4);
  GTVector<GFTYPE>          latp(4), lonp(4);
  std::vector<GFTYPE>       c0(4), sig0(4);
  std::vector<GFTYPE>       lat0(4), lon0(4); // up to 4 lumps
  std::vector<GFTYPE>       Omega;
  State            c(GDIM+1);
  GString          snut;

  PropertyTree lumpptree = ptree.getPropertyTree(sconfig);
  PropertyTree icosptree = ptree.getPropertyTree("grid_icos");
  PropertyTree advptree  = ptree.getPropertyTree("pde_burgers");
  PropertyTree nuptree   = ptree.getPropertyTree("dissipation_traits");
  GBOOL doheat   = advptree.getValue<GBOOL>("doheat");
  GBOOL bpureadv = advptree.getValue<GBOOL>("bpureadv");

  assert((bpureadv || doheat) && "Pure advection or heat must be used");

  GTVector<GTVector<GFTYPE>> *xnodes = &grid.xNodes();
  assert(grid.gtype() == GE_2DEMBEDDED && "Invalid element types");


  if ( bpureadv ) {
    for ( GSIZET j=0; j<GDIM+1; j++ ) c[j] = u[j+1];
  }

  nxy = (*xnodes)[0].size(); // same size for x, y, z

  lat0  = lumpptree.getArray<GFTYPE>("latitude0"); // lat for each lump
  lon0  = lumpptree.getArray<GFTYPE>("longitude0"); // lon for each lump
  sig0  = lumpptree.getArray<GFTYPE>("sigma"); // sig for each lump
  c0    = lumpptree.getArray<GFTYPE>("c0");  // initial concentrations for each lump
  rad   = icosptree.getValue<GFTYPE>("radius");
  E0    = lumpptree.getValue<GFTYPE>("E0");
  alpha = lumpptree.getValue<GFTYPE>("alpha",0.0);
  nlumps= lumpptree.getValue<GINT>("nlumps",1);

  alpha *= PI/180.0;

  nu    = nuptree   .getValue<GFTYPE>("nu");
  snut  = nuptree   .getValue<GString>("nu_type","constant");
  assert( snut == "constant" && "nu_type must bet set to 'constant')");

  // Convert initial locations from degrees to radians,
  // & compute initial positions of lumps in Cart coords:
  for ( GSIZET k=0; k<nlumps; k++ ) {
    lat0[k] *= PI/180.0;
    lon0[k] *= PI/180.0;
  }

  // Set velocity here. Taken to be solid body rotation,
  // u = Omega X r, where rotation rate vector, Omega
  // is computed by rotating u0 x (0, 0, 1) about x axis by
  // an amount alpha. These point to components of state u_:
  u0 = sqrt(2.0*E0);
  if ( bpureadv ) {
    for ( k=0; k<nxy; k++ ) {
      x   = (*xnodes)[0][k]; y = (*xnodes)[1][k]; z = (*xnodes)[2][k];
      r   = sqrt(x*x + y*y + z*z);
      // Compute lat & long:
      lat = asin(z/r);
      lon = atan2(y,x);
      // u_lat = u_theta = -u0 sin(lon) sin(alpha)
      // u_lon = u_phi   =  u0 (cos(theta) cos(alpha) + sin(theta)cos(lon)sin(alpha) )
      (*utmp[0])[k]  = -u0*sin(lon)*sin(alpha);
      (*utmp[1])[k]  =  u0*(cos(lat)*cos(alpha) + sin(lat)*cos(lon)*sin(alpha) );
    }
    GMTK::vsphere2cart(grid, utmp, GVECTYPE_PHYS, c);
//  GMTK::constrain2sphere(grid, c);
  }

  *u[0] = 0.0;
  for ( GSIZET k=0; k<nlumps; k++ ) {

    // Allow different sigma/concentration for each lump:
    sig  [k] = sqrt(sig0[k]*sig0[k] + 4.0*time*nu); // constant viscosity only
    si   [k] = 1.0/(sig[k]*sig[k]);
    ufact[k] = c0[k]*pow(sig0[k]/sig[k],GDIM);

    // Find where lat/lon endpoint would be at t, if alpha=0:
    latp[k]  = lat0[k];
    lonp[k]  = lon0[k] + (u0/rad) * time;

    // Find where Cart endpoint would be at t, if alpha=0:
    rt[0] = rad*cos(latp[k])*cos(lonp[k]);
    rt[1] = rad*cos(latp[k])*sin(lonp[k]);
    rt[2] = rad*sin(latp[k]);

    // Now, rotate rt about x-axis by alpha to
    // find lat/lon of final position of lump:
    xx[0] = rt[0]; xx[1] = rt[1]; xx[2] = rt[2];
    if ( time > 0 ) {
      xx[1] =  cos(alpha)*rt[1] + sin(alpha)*rt[2];
      xx[2] = -sin(alpha)*rt[1] + cos(alpha)*rt[2];
    }
    latp[k]  = asin(xx[2]/rad);
    lonp[k]  = atan2(xx[1],xx[0]);

    rt[0] = rad*cos(latp[k])*cos(lonp[k]);
    rt[1] = rad*cos(latp[k])*sin(lonp[k]);
    rt[2] = rad*sin(latp[k]);

    for ( GSIZET j=0; j<nxy; j++ ) {
      // Note: following c t is actually Integral_0^t c(t') dt', 
      //       so if c becomes a function of t, this muct change:

      x   = (*xnodes)[0][j];
      y   = (*xnodes)[1][j];
      z   = (*xnodes)[2][j];
#if 1
      r   = sqrt(x*x + y*y + z*z);
      lat = asin(z/r);
      lon = atan2(y,x);
      // Compute arclength from point to where center _should_ be:
      s     = r*acos( sin(latp[k])*sin(lat) + cos(latp[k])*cos(lat)*cos(lon-lonp[k]) );
     (*u[0])[j] += ufact[k]*exp(-s*s*si[k]);
#else
     argxp = pow(x-rt[0],2) + pow(y-rt[1],2) + pow(z-rt[2],2);
     (*u[0])[j] += ufact[k]*exp(-argxp*si[k]);
#endif
    } // end, grid-point loop

  } // end, lump loop

  return TRUE;

} // end of method impl_icosgauss




} // end, ginitstate namespace
