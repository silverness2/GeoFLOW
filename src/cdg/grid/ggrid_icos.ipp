//==================================================================================
// Description  : Object defining a (global) icosahedral grid, that
//                uses gnomonic projections to locate element vertices.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : GGrid.
//==================================================================================


//**********************************************************************************
//**********************************************************************************
// METHOD : project2sphere (1)
// DESC   : Project Cartesian coords to sphere, specified by rad argument, and
//          express in Cartesian coords.  Necessary for 2d grids.
// ARGS   : tmesh: Vector of triangles (vertices), modified to contain 
//                 projected coordinates
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::project2sphere(GTVector<GTriangle<T>> &tmesh, T rad)
{
  GString serr = "GridIcos::project2sphere (1): ";

  T          r, xlat, xlong;
  GTPoint<T> v;

  for ( GSIZET i=0; i<tmesh_.size(); i++ ) { // loop over all triangles in tmesh
    for ( GSIZET j=0; j<3; j++ ) { // guaranteed to have 3 points each
      v = *tmesh[i].v[j];
      r = v.norm();
      xlat  = asin(v.x3/r);
      xlong = atan2(v.x2,v.x1);
      xlong = xlong < 0.0 ? 2*PI+xlong : xlong;
      v.x1 = rad*cos(xlat)*cos(xlong);
      v.x2 = rad*cos(xlat)*sin(xlong);
      v.x3 = rad*sin(xlat);
      *tmesh[i].v[j] = v;
    }
  }

} // end of method project2sphere (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : project2sphere (2)
// DESC   : Project Cartesian coords to sphere, specified by rad argument.
//          Necessary for 2d grids.
// ARGS   : plist: Vector of points, modified to contain 
//                 projected coordinates. Must be 3d points.
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::project2sphere(GTVector<GTPoint<T>> &plist, T rad)
{
  GString serr = "GridIcos::project2sphere (2): ";

  T          r, xlat, xlong;
  GTPoint<T> v;

  for ( GSIZET i=0; i<plist.size(); i++ ) { // loop over all points
    v = plist[i];
    r = v.norm();
    xlat  = asin(v.x3/r);
    xlong = atan2(v.x2,v.x1);
    xlong = xlong < 0.0 ? 2*PI+xlong : xlong;
    v.x1 = rad*cos(xlat)*cos(xlong);
    v.x2 = rad*cos(xlat)*sin(xlong);
    v.x3 = rad*sin(xlat);
    plist[i] = v;
  }

} // end of method project2sphere (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : project2sphere (3)
// DESC   : Project Cartesian coords to sphere, specified by rad argument.
//          Necessary for 2d grids.
// ARGS   : plist: Vector of vectors, one vector for each dimension, modified to 
//                 contain projected coordinates. Must be 3d points.
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::project2sphere(GTVector<GTVector<T>> &plist, T rad)
{
  GString serr = "GridIcos::project2sphere (3): ";

  T r, xlat, xlong, x, y, z;

  for ( GSIZET i=0; i<plist[0].size(); i++ ) { // loop over all points
    x = plist[0][i]; y = plist[1][i]; z = plist[2][i];
    r = sqrt(x*x + y*y + z*z);
    xlat  = asin(z/r);
    xlong = atan2(y,x);
    xlong = xlong < 0.0 ? 2*PI+xlong : xlong;
    plist[0][i] = rad*cos(xlat)*cos(xlong);
    plist[1][i] = rad*cos(xlat)*sin(xlong);
    plist[2][i] = rad*sin(xlat);
  }

} // end of method project2sphere (3)


//**********************************************************************************
//**********************************************************************************
// METHOD : spherical2xyz (1)
// DESC   : Transform from spherical-polar to Cartesian coords
// ARGS   : plist: Vector of points, representing spherical coordinates 
//                 (r, lat, long), to be converted to (x, y, z), in-place
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::spherical2xyz(GTVector<GTPoint<T>*> &plist)
{
  GString serr = "GridIcos::spherical2xyz(1): ";

  T r, xlat, xlong;

  for ( GSIZET i=0; i<plist.size(); i++ ) { // loop over all points
    r = plist[i]->x1; xlat = plist[i]->x2; xlong = plist[i]->x3;
    plist[i]->x1 = r*cos(xlat)*cos(xlong);
    plist[i]->x2 = r*cos(xlat)*sin(xlong);
    plist[i]->x3 = r*sin(xlat);
  }

} // end of method spherical2xyz (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : spherical2xyz (2)
// DESC   : Transform from spherical-polar to Cartesian coords
// ARGS   : plist: Vector of points, representing spherical coordinates 
//                 (r, lat, long), to be converted to (x, y, z), in-place
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::spherical2xyz(GTVector<GTPoint<T>> &plist)
{
  GString serr = "GridIcos::spherical2xyz(2): ";

  T r, xlat, xlong;

  for ( GSIZET i=0; i<plist.size(); i++ ) { // loop over all points
    r = plist[i].x1; xlat = plist[i].x2; xlong = plist[i].x3;
    plist[i].x1 = r*cos(xlat)*cos(xlong);
    plist[i].x2 = r*cos(xlat)*sin(xlong);
    plist[i].x3 = r*sin(xlat);
  }

} // end of method spherical2xyz (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : spherical2xyz (3)
// DESC   : Transform from spherical-polar to Cartesian coords
// ARGS   : plist: Vector of points, representing spherical coordinates 
//                 (r, lat, long), to be converted to (x, y, z), in-place
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::spherical2xyz(GTVector<GTVector<T>> &plist)
{
  GString serr = "GridIcos::spherical2xyz(3): ";
  assert(plist.size() >= 3 && "Invalid dimensionality on input array");

  T r, xlat, xlong;

  for ( GSIZET i=0; i<plist[0].size(); i++ ) { // loop over all points
    r = plist[0][i]; xlat = plist[1][i]; xlong = plist[2][i];
    plist[0][i] = r*cos(xlat)*cos(xlong);
    plist[1][i] = r*cos(xlat)*sin(xlong);
    plist[2][i] = r*sin(xlat);
  }

} // end of method spherical2xyz (3)


//**********************************************************************************
//**********************************************************************************
// METHOD : xyz2spherical (1)
// DESC   : Transform from Cartesian coords to spherical-polar
// ARGS   : plist: Vector of points, representing Cartesian coordinates 
//                 (x,y,z) to be converted to (r, lat, long), in-place
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::xyz2spherical(GTVector<GTPoint<T>*> &plist)
{
  GString serr = "GridIcos::xyz2spherica(1): ";

  T r, x, y, z;

  for ( GSIZET i=0; i<plist.size(); i++ ) { // loop over all points
   x = plist[i]->x1; y = plist[i]->x2; z = plist[i]->x3;
   r = sqrt(x*x + y*y + z*z);
   plist[i]->x1 = r;
   plist[i]->x2 = asin(z/r);
   plist[i]->x3 = atan2(y,x);
   plist[i]->x3 = plist[i]->x3 < 0.0 ? 2*PI+plist[i]->x3 : plist[i]->x3;
  }

} // end of method xyz2spherical (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : xyz2spherical (2)
// DESC   : Transform from Cartesian coords to spherical-polar
// ARGS   : plist: Vector of points, representing Cartesian coordinates 
//                 (x,y,z) to be converted to (r, lat, long), in-place
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::xyz2spherical(GTVector<GTVector<T>> &plist)
{
  GString serr = "GridIcos::xyz2spherica(2): ";

  T r, x, y, z;

  for ( GSIZET i=0; i<plist.size(); i++ ) { // loop over all points
   x = plist[i][0]; y = plist[i][1]; z = plist[i][2];
   r = sqrt(x*x + y*y + z*z);
   plist[i][0] = r;
   plist[i][1] = asin(z/r);
   plist[i][2] = atan2(y,x);
   plist[i][2] = plist[i][2] < 0.0 ? 2*PI+plist[i][2] : plist[i][2];
  }

} // end of method xyz2spherical (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : xyz2spherical (3)
// DESC   : Transform from Cartesian coords to spherical-polar
// ARGS   : plist: Vector of points, representing Cartesian coordinates 
//                 (x,y,z) to be converted to (r, lat, long), in-place
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::xyz2spherical(GTVector<GTPoint<T>> &plist)
{
  GString serr = "GridIcos::xyz2spherica(3): ";

  T r, x, y, z;

  for ( GSIZET i=0; i<plist.size(); i++ ) { // loop over all points
   x = plist[i].x1; y = plist[i].x2; z = plist[i].x3;
   r = sqrt(x*x + y*y + z*z);
   plist[i].x1 = r;
   plist[i].x2 = asin(z/r);
   plist[i].x3 = atan2(y,x);
   plist[i].x3 = plist[i].x3 < 0.0 ? 2*PI+plist[i].x3 : plist[i].x3;
  }

} // end of method xyz2spherical (3)


//**********************************************************************************
//**********************************************************************************
// METHOD : cart2gnomonic
// DESC   : Transform Cartesian coords to gnomonic space
// ARGS   : clist: Vector of Cartesian points
//          rad  : radius of sphere
//          xlatc,
//          xlongc: lat and long of reference vector (usually a centroid)
//          glist: converted gnomonic points
//          
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::cart2gnomonic(GTVector<GTPoint<T>> &clist, T rad, T xlatc, T xlongc, GTVector<GTPoint<T>> &glist)
{
  GString serr = "GridIcos::cart2gnomonic: ";


  // Gnomonic transform given by (th=lat, and phi = long; thc, phic refer
  // to refernce lat and long):
  // x = rad [ cos(th) sin(phi-phic) ] / 
  //     [sin(thc)sin(th) + cos(thc)cos(th)cos(phi-phic)]
  // y = rad [ cos(thc)sin(th) - sin(thc)cos(th)cos(phi-phic) ] / 
  //     [sin(thc)sin(th) + cos(thc)cos(th)cos(th-thc)]
  // This can be rewritten as:
  // x = rad tan(longp), y = rad tan(latp) sec(longp)
  // where
  // latp  = arcsin[cos(thc)sin(th) - sin(thc)cos(th)cos(phi-phic)]
  // longp = atan [ [ cos(th) sin(phi-phic) / 
  //                  [sin(thc)sin(th) + cos(thc)cos(th)cos(phi-phic)] ] ]
  // 

  T xlatp, xlongp;
  T den, r, xlat, xlong;
  for ( GSIZET i=0; i<clist.size(); i++ ) { // loop over all points
    r      = clist[i].norm();
    xlat   = asin(clist[i].x3/r);
    xlong  = atan2(clist[i].x2,clist[i].x1);
    xlong  = xlong < 0.0 ? 2*PI+xlong : xlong;
#if 0
    den    = sin(xlatc)*sin(xlat) + cos(xlatc)*cos(xlat)*cos(xlong-xlongc);  
    xlatp  = asin( cos(xlatc)*sin(xlat) - sin(xlatc)*cos(xlat)*cos(xlong-xlongc) );
    xlongp = atan2( cos(xlat)*sin(xlong-xlongc),den );

    glist[i].x1 = rad*tan(xlongp);
    glist[i].x2 = rad*tan(xlatp)*sec(xlongp);
#else
    den    = sin(xlatc)*sin(xlat) + cos(xlatc)*cos(xlat)*cos(xlong-xlongc); 
    den    = fabs(den) < std::numeric_limits<GFTYPE>::epsilon() ? 0.0 : 1.0/den;
    glist[i].x1 = rad*cos(xlat)*sin(xlong-xlongc)*den;
    glist[i].x2 = rad*( cos(xlatc)*sin(xlat) - sin(xlatc)*cos(xlat)*cos(xlong-xlongc) ) * den;
#endif
  }

} // end of method cart2gnomonic


//**********************************************************************************
//**********************************************************************************
// METHOD : gnomonic2cart (1)
// DESC   : Transform gnomonic coords to Cartesian space
// ARGS   : glist: Vector of gnomonic coords
//          rad  : radius of sphere
//          xlatc,
//          xlongc: lat and long of reference vector (usually a centroid)
//          clist: converted Cartesian coords
//          
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::gnomonic2cart(GTVector<GTVector<T>> &glist, T rad, T xlatc, T xlongc, GTVector<GTVector<T>> &clist)
{
  GString serr = "GridIcos::gnomonic2cart (1): ";
  assert(glist.size() >= 2 && clist.size() >= 3 && "Incompaible coordinate dimensions");


  // Gnomonic transform given by (th=lat, and phi = long; thc, phic refer
  // to refernce lat and long):
  //   x = rad [ cos(th) sin(phi-phic) ] / 
  //       [sin(thc)sin(th) + cos(thc)cos(th)cos(phi-phic)]
  //   y = rad [ cos(thc)sin(th) - sin(thc)cos(th)cos(phi-phic) ] / 
  //       [sin(thc)sin(th) + cos(thc)cos(th)cos(th-thc)]
  // 
  // Reverse this here... Solving for tan(th), and cos(phi-phic), arrive at:
  //   th  = asin( cos(b)sin(thc) + y*sin(b)*cos(thc)/rho )
  //   phi = phic + atan( x*sin(b) / ( rho*cos(thc)*cos(b) - y*sin(thc)*sin(b) ) ) 
  // where
  //   rho = sqrt(x^2 + y^2)
  //   b   = atan(rho)  
  //
  // (From Wolfram Research)

  T beta, rho, x, xlat, xlong, y;
  T eps = std::numeric_limits<T>::epsilon();
  for ( GSIZET i=0; i<glist[0].size(); i++ ) { // loop over all points
    x      = glist[0][i];
    y      = glist[1][i];


    if ( fabs(x) < eps && fabs(y) < eps ) {
      xlat = xlatc;
      xlong = xlongc;
    } else {
      rho    = sqrt(x*x + y*y);
      beta   = atan(rho);
      xlat   = asin( cos(beta)*sin(xlatc) + (y*sin(beta)*cos(xlatc))/rho );
      xlong  = xlongc + atan2( x*sin(beta), 
                               rho*cos(xlatc)*cos(beta)-y*sin(xlatc)*sin(beta) );
    }

    // Convert to spherical-polar to  Cart. coordinates:
    clist[0][i] = rad*cos(xlat)*cos(xlong);
    clist[1][i] = rad*cos(xlat)*sin(xlong);
    clist[2][i] = rad*sin(xlat);
  }

} // end of method gnomonic2cart (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : gnomonic2cart (2)
// DESC   : Transform gnomonic coords to Cartesian space
// ARGS   : glist: Vector of gnomonic coords
//          rad  : radius of sphere
//          xlatc,
//          xlongc: lat and long of reference vector (usually a centroid)
//          clist: converted Cartesian coords
//          
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::gnomonic2cart(GTVector<GTPoint<T>> &glist, T rad, T xlatc, T xlongc, GTVector<GTPoint<T>> &clist)
{
  GString serr = "GridIcos::gnomonic2cart(2): ";


  // Gnomonic transform given by (th=lat, and phi = long; thc, phic refer
  // to refernce lat and long):
  //   x = rad [ cos(th) sin(phi-phic) ] / 
  //       [sin(thc)sin(th) + cos(thc)cos(th)cos(phi-phic)]
  //   y = rad [ cos(thc)sin(th) - sin(thc)cos(th)cos(phi-phic) ] / 
  //       [sin(thc)sin(th) + cos(thc)cos(th)cos(th-thc)]
  // Reverse this here... Solving for tan(th), and cos(phi-phic), arrive at:
  //   th  = asin( cos(b)sin(thc) + y*sin(b)*cos(thc)/rho )
  //   phi = phic + atan( x*sin(b) / ( rho*cos(thc)*cos(b) - y*sin(thc)*sin(b) ) ) 
  // where
  //   rho = sqrt(x^2 + y^2)
  //   b   = atan(rho)  
  // (From Wolfram Research)

  T beta, rho, x, xlat, xlong, y;
  T eps = std::numeric_limits<T>::epsilon();
  for ( GSIZET i=0; i<clist.size(); i++ ) { // loop over all points
    x      = glist[i][0];
    y      = glist[i][1];

    if ( fabs(x) < eps && fabs(y) < eps ) {
      xlat = xlatc;
      xlong = xlongc;
    } else {
      rho    = sqrt(x*x + y*y);
      beta   = atan(rho);
      xlat   = asin( cos(beta)*sin(xlatc) + (y*sin(beta)*cos(xlatc))/rho );
      xlong  = xlongc + atan2( x*sin(beta), 
                               rho*cos(xlatc)*cos(beta)-y*sin(xlatc)*sin(beta) );
    }
    // Convert to spherical-polar to  Cart. coordinates:
    clist[i][0] = rad*cos(xlat)*cos(xlong);
    clist[i][1] = rad*cos(xlat)*sin(xlong);
    clist[i][2] = rad*sin(xlat);
  }

} // end of method gnomonic2cart(2)


//**********************************************************************************
//**********************************************************************************
// METHOD : reorderverts2d
// DESC   : Reorder specified vertices to be consistent with
//          shape functions
// ARGS   : uverts : list of unordered vertices
//          overts : array of ordered vertices, returned
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::reorderverts2d(GTVector<GTPoint<T>> &uverts, GTVector<GSIZET> &isort, 
                               GTVector<GTPoint<T>> &overts)
{
  GString serr = "GridIcos::reorderverts2d: ";

  assert(uverts.size() == 4 && "Incorrect number of vertices");


  GTVector<T> x(4);
  GTVector<GSIZET> Ixy(4);


  isort.resize(4);
  for ( GSIZET i=0; i<uverts.size(); i++ ) { 
    x[i] = uverts[i].x1;
  }

  x.sortincreasing(Ixy);

  // Do 'left' -hand vertices:
  if ( uverts[Ixy[0]].x2 < uverts[Ixy[1]].x2 ) {
    isort[0] = Ixy[0];
    isort[3] = Ixy[1];
  } else {
    isort[0] = Ixy[1];
    isort[3] = Ixy[0];
  }

  // Do 'right' -hand vertices:
  if ( uverts[Ixy[2]].x2 < uverts[Ixy[3]].x2 ) {
    isort[1] = Ixy[2];
    isort[2] = Ixy[3];
  } else {
    isort[1] = Ixy[3];
    isort[2] = Ixy[2];
  }

  for ( GSIZET j=0; j<4; j++ ) overts[j] = uverts[isort[j]];
  
} // end of method reorderverts2d



//**********************************************************************************
//**********************************************************************************
// METHOD : copycast (1)
// DESC   : copy 'from' array to 'to' array, while casting. Also, make sure
//          sizes are correct.
// ARGS   : from : 'from' array
//          to   : 'to' array; may be of different type than 'from'
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::copycast(GTVector<GTVector<GFTYPE>> &from, GTVector<GTVector<T>> &to)
{

  assert(to.size() == from.size() && "Incompatible dimensions");
  for ( auto j=0; j<to.size(); j++ ) {
    to[j].resize(from[j].size()); 
    for ( auto i=0; i<to[j].size(); i++ ) {
      to[j][i] = static_cast<T>(from[j][i]);
    }
  }
  
} // end of method copycast (1)



//**********************************************************************************
//**********************************************************************************
// METHOD : copycast (2)
// DESC   : copy 'from' array to 'to' array, while casting. Also, make sure
//          sizes are correct.
// ARGS   : from : 'from' array
//          to   : 'to' array; may be of different type than 'from'
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::copycast(GTVector<GTVector<GFTYPE>*> &from, GTVector<GTVector<T>*> &to)
{

  assert(to.size() == from.size() && "Incompatible dimensions");
  for ( auto j=0; j<to.size(); j++ ) {
    to[j]->resize(from[j]->size()); 
    for ( auto i=0; i<to[j]->size(); i++ ) {
      (*to[j])[i] = static_cast<T>((*from[j])[i]);
    }
  }
  
} // end of method copycast (2)



//**********************************************************************************
//**********************************************************************************
// METHOD : copycast (3)
// DESC   : copy 'form' point to 'to' point, while casting
// ARGS   : from : 'from' point
//          to   : 'to' point; may be of different type than 'from'
// RETURNS: none.
//**********************************************************************************
template<typename T>
void GGridIcos::copycast(GTPoint<GFTYPE> &from, GTPoint<T> &to)
{

  for ( auto j=0; j<to.size(); j++ ) {
    to[j] = static_cast<T>(from[j]);
  }
  
} // end of method copycast (3)



