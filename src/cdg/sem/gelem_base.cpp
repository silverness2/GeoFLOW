//==================================================================================
// Module       : gelem_base
// Date         : 6/1/18 (DLR)
// Description  : Base class forming interfaces for all allowed 1D/2D/3D (lin/quad/hex) 
//                elements
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//               
//
//          V7 ________________ V6
//            /|     M6       /|
//       M7  / |             / |
//          /  |    F5   M5 /  |
//         /M11|  M4   F2  /   |                     z 
//     V4 /____|__________/    | M10                       y
//        |    |         | V5  |                     |     
//        |    |     M2  | F1  |                     |   /
//        | V3 |_________|_____| V2                  |  /
//    M8  | F3 / F0      |    /                      | /
//        |   /       M9 |   /                       |/________  x
//        |  / M3   F4   |  / M1                     
//        | /            | /                       
//        |/_____________|/a                     
//       V0       M0     V1
//
// Faces are labeled s.t. F0-F3 correspond to orientation of edges on bottom plane;
// F4 and F5 are bottom and top faces, respectively. Nodes are labeled implicitly starting
// from V0, increasing fastest with x, then y, then z
//
// [Note that in 2d, we use just the bottom plane.]
//==================================================================================
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "gelem_base.hpp"
#include "gshapefcn_embed.hpp"
#include "gshapefcn_hostd.hpp"
#include "gmtk.hpp"


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Default constructor
// ARGS   : none
// RETURNS: none
//**********************************************************************************
GElem_base::GElem_base()
:
bInitialized_    (FALSE),
bbasis_          (FALSE),
elemtype_        (GE_MAX),
elemid_          (0),
rootid_          (0),
nVertices_       ((GINT )pow(2.0,GDIM)),
nEdges_          (2*(GDIM-1)*GDIM),
nFaces_          (2*GDIM),
volume_          (0.0),
gshapefcn_       (NULLPTR)
{
  N_.resize(GDIM);
  gbasis_.resize(GDIM);
  gbasis_ = NULLPTR;
} // end of constructor method (1)

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (2)
// DESC   : Instantiate with basis pointers specified explicitly
// ARGS   : etype: elemet type
//          GNBasis objects for each direction
// RETURNS: none
//**********************************************************************************
GElem_base::GElem_base(GElemType etype, GNBasis<GCTYPE,GFTYPE> *b1, GNBasis<GCTYPE,GFTYPE> *b2, GNBasis<GCTYPE,GFTYPE> *b3)
:
bInitialized_    (FALSE),
elemtype_        (etype),
elemid_          (0),
rootid_          (0),
nVertices_       ((GINT )pow(2.0,GDIM)),
nEdges_          (2*(GDIM-1)*GDIM),
nFaces_          (2*GDIM),
volume_          (0.0), 
gshapefcn_       (NULLPTR)
{

  N_.resize(GDIM);
  gbasis_.resize(GDIM);
  gbasis_ = NULLPTR;

  GTVector<GNBasis<GCTYPE,GFTYPE> *> b(GDIM);

  b = NULLPTR;
  if ( b1 != NULLPTR ) b[0] = b1;
  if ( b2 != NULLPTR ) b[1] = b2;
  if ( b3 != NULLPTR ) b[2] = b3;

  if ( elemtype_ == GE_REGULAR || elemtype_ == GE_DEFORMED ) {
    gshapefcn_ = new GShapeFcn_hostd();
  }
  else if ( elemtype_ == GE_2DEMBEDDED ) { 
    gshapefcn_ = new GShapeFcn_embed();
  }

  set_basis(b);
 
} // end of constructor method (2)


//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (3)
// DESC   : Instantiate with basis pointers specified
// ARGS   : etype: elemet type
//          GNBasis object array for each direction
// RETURNS: none
//**********************************************************************************
GElem_base::GElem_base(GElemType etype, GTVector<GNBasis<GCTYPE,GFTYPE>*> &b)
:
bInitialized_    (FALSE),
elemtype_        (etype),
elemid_          (0),
rootid_          (0),
nVertices_       ((GINT )pow(2.0,GDIM)),
nEdges_          (2*(GDIM-1)*GDIM),
nFaces_          (2*GDIM),
volume_          (0.0), 
gshapefcn_       (NULLPTR)
{
  N_.resize(GDIM);
  gbasis_.resize(GDIM);
  gbasis_ = NULLPTR;

  if ( elemtype_ == GE_REGULAR || elemtype_ == GE_DEFORMED ) {
    gshapefcn_ = new GShapeFcn_hostd();
  }
  else if ( elemtype_ == GE_2DEMBEDDED ) { 
    gshapefcn_ = new GShapeFcn_embed();
  }
  set_basis(b);
 
} // end of constructor method (3)


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   : 
// ARGS   : none
// RETURNS: none
//**********************************************************************************
GElem_base::~GElem_base()
{
  if ( gshapefcn_ != NULLPTR ) delete gshapefcn_;
} // end, destructor


//**********************************************************************************
//**********************************************************************************
// METHOD : resize
// DESC   : resizes dynamically allocate quantities
//          if required
// ARGS   :
// RETURNS:  TRUE on success, else FALSE
//**********************************************************************************
void GElem_base::resize(GTVector<GINT> &newOrder)
{
  GINT  i, isame ;
  GString serr = "GElem_base::resize: ";

  assert(newOrder.size() >= GDIM && "Bad dimensionality");

  //  resize bases:
  for ( i=0, Ntot_=1; i<GDIM; i++ ) {
    gbasis_[i]->resize(newOrder[i]);
    N_     [i] = newOrder[i]+1;
    Ntot_     *= N_[i]; 
  }

  set_size();
  bInitialized_ = FALSE;

}  // end of method resize


//**********************************************************************************
//**********************************************************************************
// METHOD :  << operator method (1)
// DESC   : output stream operator
// ARGS   :
// RETURNS: ostream &
//**********************************************************************************
std::ostream &operator<<(std::ostream &str, GElem_base &e)
{
  GINT  i;
#if 0
  str << std::endl << "ElemID: " << e.elemid_ << " {" << std::endl;
  str << " ROOT ID : " << e.rootid_ ;
  str <<           " HIWORD=" << HIWORD(&(e.rootid_));
  str <<           " LOWORD=" << LOWORD(&(e.rootid_)) << std::endl;
  str << " ParentID: " << e.parentid_ ;
  str <<           " HIWORD=" << HIWORD(&(e.parentid_));
  str <<           " LOWORD=" << LOWORD(&(e.parentid_)) << std::endl;
  str << "       ID: " << e.elemid_ ;
  str <<           " HIWORD=" << HIWORD(&(e.elemid_));
  str <<           " LOWORD=" << LOWORD(&(e.elemid_)) << std::endl;
#endif
  str << " elemtype: " << e.elemtype_ ;
  str << " nVertices: " << e.nVertices_;
  str << " nEdges: " << e.nEdges_;
  str << " nFaces: " << e.nFaces_;
  str << std::endl << " lVertices: " ;
  for ( i=0; i<e.nVertices_; i++ ) str << (e.xVertices_[i]) << " "; 
  str << std::endl;
#if defined(G_IS2D)
  str << std::endl << " Edge Midpoints: " ;
  for ( i=0; i<e.nFaces_; i++ ) str << (e.spEMidpoints_[i]) << " "; 
#elif defined(G_IS3D)
  str << std::endl << " Edge Midpoints: " ;
  for ( i=0; i<e.nEdges_; i++ ) str << (e.spEMidpoints_[i]) << " "; 
  str << std::endl << " Face Midpoints: " ;
  for ( i=0; i<e.nFaces_; i++ ) str << (e.spFMidpoints_[i]) << " "; 
#endif
  str << std::endl;
  str << " bdy_indices: " << e.bdy_indices_ << std::endl;
  str << " bdy_types  : " << e.bdy_types_ << std::endl;
#if defined(G_IS3D)
//str << " face_types : " << *(e.face_types_) << std::endl;
#endif
//str << " Mask       : " <<  (e.mask_) << std::endl;
  str << std::endl << "}";

  return str;
} // end of operator <<


//**********************************************************************************
//**********************************************************************************
// METHOD : set_elemtype
// DESC   : Set element type. This must be called, together with set_basis, if
//          default constructor is used.
// ARGS   :
// RETURNS: none
//**********************************************************************************
void GElem_base::set_elemtype(GElemType etype)
{
  GString serr = "GElem_base::set_elemtype: ";
  
  elemtype_ = etype;

  if ( gshapefcn_ != NULLPTR ) delete gshapefcn_; gshapefcn_ = NULLPTR;

  if ( elemtype_ == GE_REGULAR || elemtype_ == GE_DEFORMED ) {
    gshapefcn_ = new GShapeFcn_hostd();
  }
  else if ( elemtype_ == GE_2DEMBEDDED ) {
    gshapefcn_ = new GShapeFcn_embed();
  }

} // end of method set_elemtype


//**********************************************************************************
//**********************************************************************************
// METHOD : set_basis 
// DESC   : Set basis object for coord. direction, idir in (1, 2, 3)
// ARGS   :
// RETURNS: none
//**********************************************************************************
void GElem_base::set_basis(GTVector<GNBasis<GCTYPE,GFTYPE>*> &b)
{
  GString serr = "GElem_base::set_basis: ";

  GTVector<GINT> iorder(b.size());
  for ( GINT i=0; i<GDIM; i++ ) {
    if ( b[i] == NULLPTR ) {
      std::cout << serr << "NULLPTR basis object, j=" << i << std::endl;
      exit(1);
    }
    gbasis_ [i] = b[i];
    iorder  [i] = gbasis_[i]->getOrder();
  }
  bbasis_ = TRUE;

  // resize requires shapefcn, so must set its basis first, even if it's
  // resized:
  if ( gshapefcn_ != NULLPTR ) gshapefcn_->set_basis(gbasis_);
  resize(iorder); // requires _degree_, not no. nodes; N_ will be set

} // end of method set_basis


//***********************************************************************************
//***********************************************************************************
// METHOD : set_size
// DESC   : Sizes and initializes data based on basis and GDIM of problem.
//          
// ARGS   : none.
// RETURNS: none.
//***********************************************************************************
void GElem_base::set_size()
{
  GString serr = "GElem_base::set_size: ";

  #if defined(_G_IS1D)
    set_size1d();
  #elif defined(_G_IS2D)
    set_size2d();
  #elif defined(_G_IS2D)
    set_size3d();
  #endif
  
} // end of method set_size


//***********************************************************************************
//***********************************************************************************
// METHOD : set_size1d
// DESC   : Allocates and initializes data members for use AFTER basis has been set
// ARGS   : none.
// RETURNS: none.
//***********************************************************************************
void GElem_base::set_size1d()
{
  GString serr = "GElem_base::set_size1d: ";
  assert(GDIM == 1 && "GElem_base::set_size1d: Incorrect dimensionality");

  nVertices_ = 2;
  nEdges_    = 2;
  nFaces_    = 2;

  // Vertex indices comprising edges 
  // (1 index for each edge):
  ivedge_.resize(nEdges_); 
  ivedge_[0][0] = 0;
  ivedge_[0][1] = 0;
  ivedge_[1][0] = 1;
  ivedge_[1][1] = 1;

  // Face indices defining edges
  // (no meaning in 1d or 2d)
  ifedge_.resize(0);

  // Edge indices compiring faces
  // (1 index for each face)
  ieface_.resize(0);

  // Vertex indices comprising faces
  // (1 index for each face)
  ivface_.resize(0);

  // Allocate coordinate array and vertex coordinates:
  xNodes_.resize(GDIM);
  for ( GSIZET j=0; j<GDIM; j++ ) xNodes_[j].resize(Ntot_);

  xiNodes_.resize(GDIM);
  for ( GSIZET i=0; i<GDIM; i++ )  xiNodes_[i] = gbasis_[i]->getXiNodes();;

  xVertices_.resize(nVertices_);
  
  // Allocate geometry data:
  edgeCentroid_.resize(nEdges_);
  faceCentroid_.resize(0);

#if 0
  // Metric matrix is symmetric, so we
  // use pointers so that repeated elements
  // aren't duplicated:
  dXidX_ .resize(GDIM,GDIM);
  for ( GSIZET j=0; j<GDIM; j++ ) {
    for ( GSIZET i=0; i<GDIM; i++ ) {
      dXidX_ (i,j).resize(Ntot_);
    }
  }
  Jac_.resize(Ntot_);
  faceJac_.resize(2);
  faceNormal_.resize(2);
  for ( GSIZET j=0; j<2; j++ ) {
    faceNormal_[j].resize(1);
    faceNormal_[j][0].resize(1);
  }
#endif
  mask_.resize(Ntot_);
  mask_ = 1.0;

  // Indirection indices:
  get_indirect(gbasis_, vert_indices_, edge_indices_, face_indices_);

} // end of method set_size1d


//***********************************************************************************
//***********************************************************************************
// METHOD : set_size2d
// DESC   : Initializes some quantities independent of bases, for use in construction
// ARGS   : none.
// RETURNS: none.
//***********************************************************************************
void GElem_base::set_size2d()
{
  GString serr = "GElem_base::set_size2d: ";
  assert(GDIM == 2 && "GElem_base::set_size2d: Incorrect dimensionality");

  nVertices_ = 4;
  nEdges_    = 4;
  nFaces_    = 4;


  // Vertex indices comprising edges 
  // (2 indices for each edge):
  ivedge_.resize(nEdges_); 
  for ( GSIZET j=0; j<nEdges_; j++ ) ivedge_[j].resize(2);
  ivedge_[0][0] = 0; ivedge_[0][1] = 0;
  ivedge_[1][0] = 1; ivedge_[1][1] = 2;
  ivedge_[2][0] = 3; ivedge_[2][1] = 2;
  ivedge_[3][0] = 0; ivedge_[3][1] = 3;

  // Face indices defining edges
  // (no real meaning in 2d, but we use the following)
  ifedge_.resize(nEdges_);
  for ( GSIZET j=0; j<nEdges_; j++ ) ifedge_[j].resize(1);
  ifedge_[0][0] = 0; 
  ifedge_[1][0] = 1; 
  ifedge_[2][0] = 2; 
  ifedge_[3][0] = 3; 

  // Edge indices comprising faces
  // (1 index for each face)
  ieface_.resize(nFaces_);
  for ( GSIZET j=0; j<nFaces_; j++ ) ieface_[j].resize(1);
  ieface_[0][0] = 0; 
  ieface_[1][0] = 1; 
  ieface_[2][0] = 2; 
  ieface_[3][0] = 3; 

  // Vertex indices comprising faces
  // (2 indices for each face in 2d)
  ivface_.resize(nFaces_);
  for ( GSIZET j=0; j<nFaces_; j++ ) ivface_[j].resize(2);
  ivface_[0][0] = 0; ivface_[0][1] = 1; 
  ivface_[1][0] = 1; ivface_[1][1] = 2; 
  ivface_[2][0] = 3; ivface_[2][1] = 2; 
  ivface_[3][0] = 0; ivface_[3][1] = 3; 

  // Allocate coordinate array and vertex coordinates:
  GSIZET nxy = elemtype_ == GE_2DEMBEDDED ? GDIM+1: GDIM;
  xNodes_.resize(nxy);
  for ( GSIZET j=0; j<nxy; j++ ) xNodes_[j].resize(Ntot_);

  xiNodes_.resize(GDIM);
  for ( GSIZET i=0; i<GDIM; i++ )  xiNodes_[i] = gbasis_[i]->getXiNodes();;

#if 0
  xiNodes_.resize(Ntot_);
  GSIZET n;
  GTVector<GFTYPE> *xinodes1d[GDIM];
  for ( GSIZET i=0; i<GDIM; i++ )  xinodes1d = gbasis_[i];
  for ( GSIZET j=0, n=0; j<N_[0]; j++ ) {
    for ( GSIZET i=0; i<N_[0]; i++, n++ ) {
      xiNodes_[n] = (*xinodes1d[0])[i] * (*xinodes1d[1])[j];
    }
  }
#endif

  xVertices_.resize(nVertices_);
  
  // Allocate geometry data:
  edgeCentroid_.resize(nEdges_);
  faceCentroid_.resize(nEdges_);
  nxy = elemtype_ == GE_2DEMBEDDED ? GDIM+1: GDIM;

#if 0
  // Metric matrix is symmetric, so we
  // use pointers so that repeated elements
  // aren't duplicated:
  dXidX_ .resize(nxy,nxy);
  for ( GSIZET j=0; j<nxy; j++ ) {
    for ( GSIZET i=0; i<nxy; i++ ) {
      dXidX_ (i,j).resize(Ntot_);
    }
  }
  Jac_.resize(Ntot_);

  faceJac_.resize(nFaces_);
  faceJac_[0].resize(N_[0]);
  faceJac_[1].resize(N_[1]);
  faceJac_[2].resize(N_[0]);
  faceJac_[3].resize(N_[1]);

  faceNormal_.resize(nFaces_);
  if ( elemtype_ == GE_2DEMBEDDED ) {
    for ( GSIZET j=0; j<nFaces_; j++ ) {
      faceNormal_[j].resize(3);
      for ( GSIZET k=0; k<3; k++ ) faceNormal_[j][k].resize(N_[j%2]);
    }
  }
  else {
    for ( GSIZET j=0; j<nFaces_; j++ ) {
      faceNormal_[j].resize(2);
      for ( GSIZET k=0; k<2; k++ ) faceNormal_[j][k].resize(N_[j%2]);
    }
  }
#endif
 
  mask_.resize(Ntot_);
  mask_ = 1.0;
  
  // Indirection indices:
  get_indirect(gbasis_, vert_indices_, edge_indices_, face_indices_);


} // end of method set_size2d


//***********************************************************************************
//***********************************************************************************
// METHOD : set_size3d
// DESC   : Initializes some quantities independent of bases, for use in construction
// ARGS   : none.
// RETURNS: none.
//***********************************************************************************
void GElem_base::set_size3d()
{
  GString serr = "GElem_base::set_size3d: ";
  assert(GDIM == 3 && "GElem_base::set_size3d: Incorrect dimensionality");

  GINT   i, j;

  nVertices_ = 8;
  nEdges_    = 12;
  nFaces_    = 6;


  // Vertex indices comprising edges 
  ivedge_.resize(nEdges_); 
  for ( GSIZET j=0; j<nEdges_; j++ ) ivedge_[j].resize(2);
  ivedge_ [0][0] = 0; ivedge_ [0][1] = 1; 
  ivedge_ [1][0] = 1; ivedge_ [1][1] = 2; 
  ivedge_ [2][0] = 3; ivedge_ [2][1] = 2; 
  ivedge_ [3][0] = 0; ivedge_ [3][1] = 3; 
  ivedge_ [4][0] = 4; ivedge_ [4][1] = 5; 
  ivedge_ [5][0] = 5; ivedge_ [5][1] = 6; 
  ivedge_ [6][0] = 7; ivedge_ [6][1] = 6; 
  ivedge_ [7][0] = 4; ivedge_ [7][1] = 7; 
  ivedge_ [8][0] = 0; ivedge_ [8][1] = 4; 
  ivedge_ [9][0] = 1; ivedge_ [9][1] = 5; 
  ivedge_[10][0] = 2; ivedge_[10][1] = 6; 
  ivedge_[11][0] = 3; ivedge_[11][1] = 7; 

  // Face indices defining edges
  ifedge_.resize(nEdges_);
  for ( GSIZET j=0; j<nEdges_; j++ ) ivedge_[j].resize(1);
  ifedge_ [0][0] = 0; ifedge_ [0][1] = 4; 
  ifedge_ [1][0] = 1; ifedge_ [1][1] = 4; 
  ifedge_ [2][0] = 2; ifedge_ [2][1] = 4; 
  ifedge_ [3][0] = 3; ifedge_ [3][1] = 4; 
  ifedge_ [4][0] = 0; ifedge_ [4][1] = 5; 
  ifedge_ [5][0] = 1; ifedge_ [5][1] = 5; 
  ifedge_ [6][0] = 2; ifedge_ [6][1] = 5; 
  ifedge_ [7][0] = 3; ifedge_ [7][1] = 5; 
  ifedge_ [8][0] = 3; ifedge_ [8][1] = 0; 
  ifedge_ [9][0] = 0; ifedge_ [9][1] = 1; 
  ifedge_[10][0] = 1; ifedge_[10][1] = 2; 
  ifedge_[11][0] = 2; ifedge_[11][1] = 3; 

  // Edge indices comprising faces
  ieface_.resize(nFaces_);
  for ( GSIZET j=0; j<nFaces_; j++ ) ieface_[j].resize(4);
  ieface_[0][0] =  0; ieface_[0][1] =  9; ieface_[0][2] =  4; ieface_[0][3] =  8;
  ieface_[1][0] =  1; ieface_[1][1] = 10; ieface_[1][2] =  5; ieface_[1][3] =  9;
  ieface_[2][0] =  2; ieface_[2][1] = 10; ieface_[2][2] =  6; ieface_[2][3] = 11;
  ieface_[3][0] =  3; ieface_[3][1] = 11; ieface_[3][2] =  7; ieface_[3][3] =  8;
  ieface_[4][0] =  0; ieface_[4][1] =  1; ieface_[4][2] =  2; ieface_[4][3] =  3;
  ieface_[5][0] =  4; ieface_[5][1] =  5; ieface_[5][2] =  6; ieface_[5][3] =  7;

  // Vertex indices comprising faces
  ivface_.resize(nFaces_);
  for ( GSIZET j=0; j<nFaces_; j++ ) ivface_[j].resize(4);
  ivface_[0][0] = 0; ivface_[0][1] = 1; ivface_[0][2] = 5; ivface_[0][3] = 4;
  ivface_[1][0] = 1; ivface_[1][1] = 2; ivface_[1][2] = 6; ivface_[1][3] = 5;
  ivface_[2][0] = 3; ivface_[2][1] = 2; ivface_[2][2] = 6; ivface_[2][3] = 7;
  ivface_[3][0] = 0; ivface_[3][1] = 3; ivface_[3][2] = 7; ivface_[3][3] = 4;
  ivface_[4][0] = 0; ivface_[4][1] = 1; ivface_[4][2] = 2; ivface_[4][3] = 3;
  ivface_[5][0] = 4; ivface_[5][1] = 5; ivface_[5][2] = 6; ivface_[5][3] = 7;


  // Allocate coordinate array and vertex coordinates:
  xNodes_.resize(GDIM);
  for ( GSIZET j=0; j<GDIM; j++ ) xNodes_[j].resize(Ntot_);

  xiNodes_.resize(GDIM);
  for ( GSIZET i=0; i<GDIM; i++ )  xiNodes_[i] = gbasis_[i]->getXiNodes();;

#if 0
  xiNodes_.resize(Ntot_);
  GSIZET n;
  GTVector<GFTYPE> *xinodes1d[GDIM];
  for ( GSIZET i=0; i<GDIM; i++ )  xinodes1d = gbasis_[i];
   for ( GSIZET k=0, n=0; n<N_[2]; k++ ) {
    for ( GSIZET j=0; j<N_[0]; j++ ) {
      for ( GSIZET i=0; i<N_[0]; i++, n++ ) {
        xiNodes_[n] = (*xinodes1d[0])[i] * (*xinodes1d[1])[i] * (*xinodes1d[2])[k];
      }
    }
#endif

  xVertices_.resize(nVertices_);
  
  // Allocate geometry data:
  edgeCentroid_.resize(nEdges_);
  faceCentroid_.resize(nEdges_);

#if 0
  // Metric matrix is symmetric, so we
  // use pointers so that repeated elements
  // aren't duplicated:
  dXidX_ .resize(GDIM,GDIM);
  for ( GSIZET j=0; j<GDIM; j++ ) {
    for ( GSIZET i=0; i<GDIM; i++ ) {
      dXidX_ (i,j).resize(Ntot_);
    }
  }
  Jac_.resize(Ntot_);

  faceJac_.resize(nFaces_);
  faceJac_[0].resize(N_[0]*N_[2]);
  faceJac_[1].resize(N_[1]*N_[2]);
  faceJac_[2].resize(N_[0]*N_[2]);
  faceJac_[3].resize(N_[1]*N_[2]);
  faceJac_[4].resize(N_[0]*N_[1]);
  faceJac_[5].resize(N_[0]*N_[1]);

  faceNormal_.resize(nFaces_);
  for ( GSIZET j=0; j<nFaces_; j++ ) {
    faceNormal_[j].resize(3);
    if ( j < 4 ) {
      for ( GSIZET k=0; k<3; k++ ) faceNormal_[j][j].resize(N_[j%2]*N_[2]);
    }
    else {
      for ( GSIZET k=0; k<3; k++ ) faceNormal_[j][j].resize(N_[0]*N_[1]);
    }
  }
#endif

  mask_.resize(Ntot_);
  mask_ = 1.0;
  
  // Indirection indices:
  get_indirect(gbasis_, vert_indices_, edge_indices_, face_indices_);

} // end of method set_size3d


//***********************************************************************************
//***********************************************************************************
// METHOD : build_elem
// DESC   : Compute all geometry, and other related quantities. Call must have
//          been made to set_size prior to entry.
// ARGS   : none.
// RETURNS: none.
//***********************************************************************************
void GElem_base::build_elem()
{
  GString serr = "GElem_base::build_elem: ";

#if defined(_G_IS1D )
  build_elem1d();
#elif defined (_G_IS2D )
  build_elem2d();
#elif defined (_G_IS3D )
  build_elem3d();
#endif

} // end of method build_elem


//***********************************************************************************
//***********************************************************************************
// METHOD : build_elem1d
// DESC   : Compute all geometry, and other related quantities. Assumes
//          xNodes have been set and that set_size has been called (after
//          basis is set)
// ARGS   : none.
// RETURNS: none.
//***********************************************************************************
void GElem_base::build_elem1d()
{
  GString serr = "GElem_base::build_elem1d: ";
  assert(GDIM == 1 && "GElem_base::build_elem1d: Incorrect dimensionality");

  GSIZET i, j, k;
  GFTYPE xi1p, xi1m;
  GFTYPE xi2p, xi2m;

  // With grid x values ('xNodes') set, find vertices:
  for ( j=0; j<nVertices_; j++ ) {
    for ( k=0; k<xNodes_.size(); k++ ) xVertices_[j][k] = xNodes_[k][vert_indices_[j][0]];
  }

#if 0
  // Compute metrics, Jacobians, etc:
  dogeom1d(gbasis_, dXidX_, Jac_, faceJac_);
#endif

  // Compute edge/face centroids:
  // NOTE: not currently used
  edgeCentroid_[0] = xVertices_[0] ;
  edgeCentroid_[1] = xVertices_[1] ;
  for ( k=0; k<nFaces_; k++ ) faceCentroid_[k] = edgeCentroid_[k];
  
  // Compute element centroid:
  elemCentroid_ = (xVertices_[0] + xVertices_[1])  * 0.5;
  

} // end of method build_elem1d


//***********************************************************************************
//***********************************************************************************
// METHOD : build_elem2d
// DESC   : Compute all geometry, and other related quantities. Assumes
//          xNodes have been set and that set_size has been called (after
//          basis is set)
// ARGS   : none.
// RETURNS: none.
//***********************************************************************************
void GElem_base::build_elem2d()
{
  GString serr = "GElem_base::build_elem2d: ";
  assert(GDIM == 2 && "GElem_base::build_elem2d: Incorrect dimensionality");

  GSIZET i, j, k, n;

  // With grid x values ('xNodes') set, find vertices:
  GTPoint<GFTYPE> pt;
  for ( j=0; j<nVertices_; j++ ) {
    for ( k=0; k<xNodes_.size(); k++ ) {
      pt[k]  = xNodes_[k][vert_indices_[j][0]];
      xVertices_[j][k] = pt[k];
    }
  }

#if 0
  // Compute metrics, Jacobians, etc:
  dogeom2d(gbasis_, dXidX_, Jac_, faceJac_, edge_indices_);
#endif
  
  // Compute edge/face centroids:
  edgeCentroid_[0] = (xVertices_[0] + xVertices_[1]) * 0.5;
  edgeCentroid_[1] = (xVertices_[1] + xVertices_[2]) * 0.5;
  edgeCentroid_[2] = (xVertices_[3] + xVertices_[2]) * 0.5;
  edgeCentroid_[3] = (xVertices_[0] + xVertices_[3]) * 0.5;
  for ( k=0; k<nFaces_; k++ ) faceCentroid_[k] = edgeCentroid_[k];
  
  // Compute element centroid:
  elemCentroid_ = (xVertices_[0] + xVertices_[1] + xVertices_[2] + xVertices_[3]) * 0.25;
  

} // end of method build_elem2d


//***********************************************************************************
//***********************************************************************************
// METHOD : build_elem3d
// DESC   : Compute all geometry, and other related quantities. Assumes
//          xNodes have been set and that set_size has been called (after
//          basis is set)
// ARGS   : none.
// RETURNS: none.
//***********************************************************************************
void GElem_base::build_elem3d()
{
  GString serr = "GElem_base::build_elem3d: ";
  assert(GDIM == 3 && "GElem_base::build_elem3d: Incorrect dimensionality");

  GSIZET i, j, k, n;

  // With grid x values ('xNodes') set, find vertices:
  for ( j=0; j<nVertices_; j++ ) {
    for ( k=0; k<xNodes_.size(); k++ ) xVertices_[j][k] = xNodes_[k][vert_indices_[j][0]];
  }

#if 0
  // Compute metrics, Jacobians, etc:
  dogeom3d(gbasis_, dXidX_, Jac_, faceJac_, face_indices_);
#endif

  // Compute edge/face centroids:
  // NOTE: not currently used
  for ( j=0; j<4; j++ ) { // bottom plane
    edgeCentroid_ [j] = (xVertices_[j] + xVertices_[(j+1)%4]) * 0.5;
  }
  for ( j=0; j<4; j++ ) { // top plane
    edgeCentroid_ [j+4] = (xVertices_[j+4] + xVertices_[(j+1)%4+4]) * 0.5;
  }
  for ( j=0; j<4; j++ ) { // vertical edges
    edgeCentroid_ [j+8] = (xVertices_[j] + xVertices_[(j+1)%4+4]) * 0.5;
  }
  faceCentroid_[0] = (xVertices_[0] + xVertices_[1] + xVertices_[4] + xVertices_[5]) * 0.25;
  faceCentroid_[1] = (xVertices_[1] + xVertices_[2] + xVertices_[5] + xVertices_[6]) * 0.25;
  faceCentroid_[2] = (xVertices_[2] + xVertices_[3] + xVertices_[6] + xVertices_[7]) * 0.25;
  faceCentroid_[3] = (xVertices_[0] + xVertices_[3] + xVertices_[4] + xVertices_[7]) * 0.25;
  faceCentroid_[4] = (xVertices_[0] + xVertices_[1] + xVertices_[2] + xVertices_[3]) * 0.25;
  faceCentroid_[5] = (xVertices_[4] + xVertices_[5] + xVertices_[6] + xVertices_[7]) * 0.25;
  
  // Compute element centroids:
  elemCentroid_ = 0.0;
  for ( j=0; j<nVertices_; j++ ) { // average vertices
    elemCentroid_ += xVertices_[j] ;
  }
  elemCentroid_ *= 1.0/static_cast<GFTYPE>(nVertices_);
  

} // end of method build_elem3d


//***********************************************************************************
//***********************************************************************************
// METHOD : dogeom1d
// DESC   : Compute all geometry objects that depend on basis in 1d
//          NOTE: no resizing here, as all return quantities may be
//                global here, and simply restricted to its range
//                of this element.
// ARGS   : 
//          rij  : dx^j/dxi^i, computed here, but otherwise is temp space
//                 that may be reallocated
//          irij : dxi^j/dx^i matrix to be created; may be global
//          jac  : Jacobian to be created ; gmay be lobal
//          fjac : face Jacobians to be created; may be global 
//          faceNormal  : normal at face nodes; may be global 
// RETURNS: none.
//***********************************************************************************
void GElem_base::dogeom1d(GTMatrix<GTVector<GFTYPE>> &rij, GTMatrix<GTVector<GFTYPE>> &irij, GTVector<GFTYPE> &jac, GTVector<GFTYPE> &fjac, GTVector<GTVector<GFTYPE>> &faceNormal)
{
  GString serr = "GElem_base::dogeom1d: ";
  assert(gshapefcn_ != NULLPTR && "GElem_base::dogeom1d: No shape function specified");

  GSIZET i, j, k, l, m, nnodes;
  GTVector<GTVector<GFTYPE>*>  xi_ev(gbasis_.size()); // ref points at which to evaluate shape fcns

  // Get total no. nodes to evaluate at:
  for ( k=0,nnodes=1; k<gbasis_.size(); k++ ) nnodes *= (gbasis_[k]->getOrder()+1);
  assert(nnodes != 1 && "GElem_base::dogeom1d: Evaluation basis not set");

  // Get reference node points:
  for ( j=0; j<gbasis_.size(); j++ ) xi_ev[j] = gbasis_[j]->getXiNodes();


  // Compute metric elements: Rij = dx^i/dxi^j, and Jacobian:
  // Spatial points are given in high order as:
  //       x^j = Sum_k x^j_k N_k(xi,eta, ...),
  // where N_k is the shape function. To compute Rij, we
  // just need the derivative of N_i, wrt each xi-coordinate:
  GTVector<GFTYPE> dNi(Ntot_);// shape function derivative
  GTVector<GINT>   I(1);      // tensor product 'index'

  rij .resizem(GDIM,GDIM);
  for ( m=0; m<GDIM; m++ ) {
    for ( l=0; l<GDIM; l++ ) {
      rij(l,m) = 0.0;
      for ( i=0; i<N_[0]; i++ ) {
        I[0] = i;
        gshapefcn_->dNdXi(I, m+1, xi_ev, dNi);
        rij(l,m) += dNi*xNodes_[0][i];
      } // i-loop
    } // l-loop
  } // m-loop
  jac = rij(0,0);

  for ( j=0; j<Ntot_; j++ ) {
    irij(0,0)[j] = 1.0 / rij(0.0)[j];
  }

#if 0
  // Compute metric: gij = Sum_k dxi^i/dx^k dxi^j/dx^k * Jac * W:
  GTVector<GTVector<GFTYPE>*> W(GDIM);
  for ( GSIZET j=0; j<GDIM; j++ ) {
    W[i]= gbasis_[j]->getWeights();
  }
  for ( j=0; j<Ntot_; j++ ) {
    gij(0,0)[j] = (*W[0])[j] / rij(0.0)[j];
  }

  // Compute face Jacobians:
  fjac[0][0] = rij(0,0)[0];
  fjac[1][0] = rij(0,0)[Ntot_-1];
#endif

  // Compute edge/face normals: 
  // A vector of GFPoints for each face:
  faceNormal[0][0] = -1.0;
  faceNormal[0][1] =  1.0;

} // end of method dogeom1d


//***********************************************************************************
//***********************************************************************************
// METHOD : dogeom2d
// DESC   : Compute all geometry objects that depend on basis in 2d
//          NOTE: no resizing here, as all return quantities may be
//                global here, and simply restricted to its range
//                of this element.
// ARGS   : 
//          rij  : dx^j/dxi^i, computed here, but otherwise is temp space
//                 that may be reallocated
//          irij : dxi^j/dx^i matrix to be created; may be global
//          jac  : Jacobian to be created ; gmay be lobal
//          fjac : face Jacobians to be created; may be global 
//          faceNormal  : normal at face nodes; may be global 
// RETURNS: none.
//***********************************************************************************
void GElem_base::dogeom2d(GTMatrix<GTVector<GFTYPE>> &rij, GTMatrix<GTVector<GFTYPE>> &irij, GTVector<GFTYPE> &jac, GTVector<GFTYPE> &fjac, GTVector<GTVector<GFTYPE>> &faceNormal)
{
  GString serr = "GElem_base::dogeom2d: ";
  assert(gshapefcn_ != NULLPTR && "GElem_base::dogeom2d: No shape function specified");

  GSIZET i, j, k, l, m, n, nnodes;
  GTVector<GINT> N(2);
  GTVector<GFTYPE> L(2);
  GTVector<GTVector<GFTYPE>*>  xi_ev(gbasis_.size()); // ref points at which to evaluate shape fcns

  for ( k=0,nnodes=1; k<gbasis_.size(); k++ ) {
    N[k] = gbasis_[k]->getOrder()+1;
    nnodes *= N[k];
  }
  assert(nnodes > 1 && "GElem_base::dogeom2d: Evaluation basis not set");

  // Set reference node points:
  for ( j=0; j<gbasis_.size(); j++ ) xi_ev[j] = gbasis_[j]->getXiNodes();

  // Get box length from vertices:
  L[0] = xVertices_[1].x1 - xVertices_[0].x1;
  L[1] = xVertices_[2].x2 - xVertices_[1].x2;

  // Compute metric elements: Rij = dx^i/dxi^j, and Jacobian.
  // Spatial points are given in high order as:
  //       x^j = Sum_K x^j_k N_K(xi,eta, ...),
  // where N_K is the shape function. To compute Rij, we
  // just need the derivative of N_I, wrt each xi-coordinate:
  GBOOL   pChk;        // check for positive-difiniteness
  GTVector<GFTYPE> dNi(nnodes);// shape function derivative
  GTVector<GFTYPE> tmp(nnodes);// tmp space
  GTVector<GINT>   I(2);      // tensor produc index

  // Can have 'embedded' coords, so # Cartesian coordinates may be > GDIM;
  // but the total number of node points in each metrix element will 
  // still be (h1-order+1) X (h2-order+1):
  GSIZET nxy = elemtype_ == GE_2DEMBEDDED ? GDIM+1: GDIM;
  if ( elemtype_ == GE_2DEMBEDDED ) {
    rij.resize(nxy,nxy);
//  irij.resize(nxy,nxy);
    for ( l=0; l<nxy; l++ ) { // rij matrix element col
      for ( k=0; k<nxy; k++ ) { // rij matrix element row
        rij(k,l).resize(Ntot_);
//      irij(k,l).resize(Ntot_);
        tmp  = 0.0;
        for ( j=0, n=0; j<gbasis_[1]->getOrder()+1; j++ ) { // evaluate gbasis at xi_ev
          for ( i=0; i<gbasis_[0]->getOrder()+1; i++, n++ ) { // evaluate gbasis at xi_ev
            I[0] = i; I[1] = j;
            gshapefcn_->dNdXi(I, l+1, xi_ev, dNi); // xi^m-th deriv of shape function I
            tmp += dNi*xNodes_[k][n];  // multiply by spatial coord
          } // i-loop
        } // j-loop
        rij(k,l) = tmp;
      } // k-loop
    } // l-loop
  } else {  // dXi/dX are just constants for reg elements:
//  irij.resize(nxy,1);
    for ( k=0; k<nxy; k++ ) { // rij matrix element col
//    irij(k,0).bconstdata(TRUE);
      irij(k,0) = 2.0/L[k];
    } // k-loop
  }


  // Compute Jacobian; test for positive-definiteness:
#if 0
  if ( elemtype_ == GE_2DEMBEDDED ) {
    Jac_embed(rij, jac, pChk, NULLPTR, 0);
  }
  else {
    det (rij, jac, pChk, NULLPTR, 0);
  }
#else
  if ( elemtype_ == GE_2DEMBEDDED ) {
    det (rij, jac, pChk, NULLPTR, 0);
  }
  else if ( elemtype_ == GE_REGULAR ) {
//  jac.bconstdata(TRUE);
    jac[0] = 0.25*L[0]*L[1];
  }
#endif

  // Find inverse of rij:
  inv(rij, jac, irij);

#if 0
  // gij moved to GHelmholtz operator, only place (so far) where
  // they are required.
  GTVector<GTVector<GFTYPE>*> W(GDIM);
  for ( GSIZET j=0; j<GDIM; j++ ) {
    W[j]= gbasis_[j]->getWeights();
  }

  // Create metric elements:
  //   gij = Sum_k dxi^i/dx^k dxi^j/dx^k * Jac * W:
  for ( j=0; j<nxy; j++ ) { // gij matrix element col
    for ( i=0; i<nxy; i++ ) { // gij matrix element row
       for ( k=0; k<nxy; k++ ) { 
         for ( m=0, n=0; m<N[1]; m++ ) {
           for ( l=0; l<N[0]; l++,n++ ) {
             (*gij(i,j))[n] = irij(i,k)[n] * irij(j,k)[n]
                            * (*W[0])[l] * (*W[1])[m] * jac[n];
           }
         }
       }
    }
  }
#endif

//assert(pChk && "Jacobian not positive definite");

  // Compute face Jacobians. Linearize edge nodes, including them
  // in order:
  GINT ntot=0;
  GTVector<GINT> iedge;
  for ( j=0; j<nEdges_; j++ ) ntot += edge_indices_[j].size();
  iedge.resize(ntot);
  for ( j=0,ntot=0; j<nEdges_; j++ ) { 
    for ( k=0; k<edge_indices_[j].size(); k++ ) 
      iedge[ntot++] = edge_indices_[j][k];
  }
  det(rij, fjac, pChk, iedge.data(), iedge.size()); 

  // Compute edge/face normals: 
  set_faceNormal2d(rij, faceNormal);

} // end of method dogeom2d


//***********************************************************************************
//***********************************************************************************
// METHOD : dogeom3d
// DESC   : Compute all geometry objects that depend on basis in 3d
//          NOTE: no resizing here, as all return quantities may be
//                global here, and simply restricted to its range
//                of this element.
// ARGS   : 
//          rij  : dx^j/dxi^i, computed here, but otherwise is temp space
//                 that may be reallocated
//          irij : dxi^j/dx^i matrix to be created; may be global
//          jac  : Jacobian to be created ; gmay be lobal
//          fjac : face Jacobians to be created; may be global 
//          faceNormal  : normal at face nodes; may be global 
// RETURNS: none.
//***********************************************************************************
void GElem_base::dogeom3d(GTMatrix<GTVector<GFTYPE>> &rij, GTMatrix<GTVector<GFTYPE>> &irij, GTVector<GFTYPE> &jac, GTVector<GFTYPE> &fjac, GTVector<GTVector<GFTYPE>> &faceNormal)
{
  GString serr = "GElem_base::dogeom3d: ";
  assert(gshapefcn_ != NULLPTR && "GElem_base::dogeom3d: No shape function specified");

  GSIZET i, j, k, l, m, n, p, nnodes;
  GTVector<GINT> N(2);
  GTVector<GFTYPE> L(3);
  GTVector<GTVector<GFTYPE>*>  xi_ev(gbasis_.size()); // ref points at which to evaluate shape fcns

  for ( k=0,nnodes=1; k<gbasis_.size(); k++ ) {
    N[k] = gbasis_[k]->getOrder()+1;
    nnodes *= N[k];
  }

  assert(nnodes > 1 && "GElem_base::dogeom3d: Evaluation basis not set");

  // Get box length from vertices:
  L[0] = xVertices_[1].x1 - xVertices_[0].x1;
  L[1] = xVertices_[2].x2 - xVertices_[1].x2;
  L[2] = xVertices_[4].x3 - xVertices_[0].x3;

  // Set reference node points:
  for ( j=0; j<gbasis_.size(); j++ ) xi_ev[j] = gbasis_[j]->getXiNodes();


  // Compute metric elements: Rij = dx^i/dxi^j, and Jacobian:
  // Spatial points are given in high order as:
  //       x^j = Sum_k x^j_k N_K(xi,eta, ...),
  // where N_K is the shape function. To compute Rij, we
  // just need the derivative of N_I, wrt each xi-coordinate:
  GBOOL   pChk;        // check for positive-difiniteness
  GTVector<GFTYPE> dNi(nnodes);// shape function derivative
  GTVector<GINT>   I(2);      // tensor produc index

  // Can have 'embedded' coords, so # Cartesian coordinates may not be GDIM;
  // but the total number of node points in each metrix element will 
  // still be (h1-order+1) X (h2-order+1):
  GSIZET nxy = elemtype_ == GDIM;
  if ( elemtype_ == GE_DEFORMED ) {
    rij.resizem(nxy,nxy);
//  irij.resize(nxy,nxy);
    for ( m=0; m<nxy; m++ ) { // matrix element col
      for ( l=0; l<nxy; l++ ) { // matrix element row
        rij(k,l).resizem(Ntot_);
//      irij(k,l).resize(Ntot_);
        rij(l,m) = 0.0;
        for ( k=0, n=0; k<gbasis_[2]->getOrder()+1; k++ ) {
          for ( j=0; j<gbasis_[1]->getOrder()+1; j++ ) {
            for ( i=0; i<gbasis_[0]->getOrder()+1; i++, n++ ) {
              I[0] = i; I[1] = j; I[2] = k;
              gshapefcn_->dNdXi(I, m+1, xi_ev, dNi); // m-th deriv of shape function I
              rij(l,m) += dNi*xNodes_[l][n]; // multiply by spatial coord
            } // i-loop
          } // j-loop
        } // k-loop
      } // l-loop
    } // m-loop
  } else {  // dXi/dX are just constants for reg elements:
//  irij.resize(nxy,1);
    for ( k=0; k<nxy; k++ ) { 
//    irij(k,0).bconstdata(TRUE);
      irij(k,0) = 2.0/L[k];
    } // k-loop
  }

  // Compute Jacobian; test for positive-definiteness:
  if ( elemtype_ == GE_DEFORMED ) {
    det (rij, jac, pChk, NULLPTR, 0);
  }
  else if ( elemtype_ == GE_REGULAR ) {
//  jac.bconstdata(TRUE);
    jac[0] = 0.125*L[0]*L[1]*L[2];
  }

  assert(pChk && "Jacobian not positive definite");

  inv(rij, jac, irij);

#if 0
  GTVector<GTVector<GFTYPE>*> W(GDIM);
  for ( GSIZET j=0; j<GDIM; j++ ) {
    W[i]= gbasis_[j]->getWeights();
  }
  // Create metric elements:
  //   gij = Sum_k dxi^i/dx^k dxi^j/dx^k * Jac * W:
  for ( j=0; j<nxy; j++ ) { // gij matrix element col
    for ( i=0; i<nxy; i++ ) { // gij matrix element row
      for ( k=0; k<nxy; k++ ) { 
        for ( p=0, n=0; p<N[2]; p++ ) {
          for ( m=0; m<N[1]; m++ ) {
            for ( l=0; l<N[0]; l++,n++ ) {
              (*gij(i,j))[n] = irij(i,k)[n] * irij(j,k)[n]
                             * (*W[0])[l] * (*W[1])[m] * (*W[2])[p] * jac[n];
            }
          }
        }
      }
    }
  }
#endif

  // Compute face Jacobians. Linearize edge nodes, including them
  // in order:
  GINT ntot=0;
  GTVector<GINT> iface;
  for ( j=0; j<nFaces_; j++ ) ntot += face_indices_[j].size();
  iface.resize(ntot);
  for ( j=0,ntot=0; j<nFaces_; j++ ) { 
    for ( k=0; k<face_indices_[j].size(); k++ ) 
      iface[ntot++] = face_indices_[j][k];
  }
  det(rij, fjac, pChk, iface.data(), iface.size()); 

  // Compute face normals: 
  set_faceNormal3d(rij, faceNormal);

} // end of method dogeom3d


//***********************************************************************************
//***********************************************************************************
// METHOD : det
// DESC   : Compute in vectorized way the determinant of specified 
//          (vectorized) matrix
// ARGS   : 
//          G    : matrix to find determinant of. Each element of G is assumed
//                 to have the same length. 
//          detv  : determinant computed. Must be at least the same length of
//                 each element of G, of size nind, if pind is non-NULLPTR
//          pChk : are all elements positive definite?
//          pind : indirection indices, used if non-NULLPTR
//          nind : number of indirection indices in pind
// RETURNS: none.
//***********************************************************************************
void GElem_base::det(GMVFType &G, GTVector<GFTYPE> &detv, GBOOL &pChk, GINT *pind, GINT nind )
{
  GString serr = "GElem_base::det: ";

  assert((G.size(1) == 2 && G.size(2) == 2) || 
         (G.size(1) == 3 || G.size(2) == 3) && "Invalid matrix dimension");

  // Compute det, and check Jacobian for being positive-definiteness:
  GSIZET n, k;
  pChk = TRUE;
  if ( G.size(1) == 2 ) { // 2x2 matrix:
    if ( pind != NULLPTR ) {
      for ( k=0; k<nind; k++ ) { // loop over desired indices only
        n = pind[k];
        detv[k] = G(0,0)[n]*G(1,1)[n] - G(0,1)[n]*G(1,0)[n];
        pChk = pChk && fabs(detv[k]) > fabs(std::numeric_limits<GFTYPE>::epsilon());  // test for zero Jac 
      }
    }
    else  {
      for ( n=0; n<detv.size(); n++ ) {
        detv[n] = G(0,0)[n]*G(1,1)[n] - G(0,1)[n]*G(1,0)[n];
        pChk = pChk && fabs(detv[n]) > fabs(std::numeric_limits<GFTYPE>::epsilon());  // test for zero Jac 
      }
    }
    return;
  }
  
  if ( pind != NULLPTR ) { // 3x3 matrix:
    for ( k=0; k<nind; k++ ) { // loop over desired indices only
      n = pind[k];
      detv[k] = G(0,0)[n]*(G(1,1)[n]*G(2,2)[n] - G(1,2)[n]*G(2,1)[n])
              - G(0,1)[n]*(G(1,0)[n]*G(2,2)[n] - G(1,2)[n]*G(2,0)[n])
              + G(0,2)[n]*(G(1,0)[n]*G(2,1)[n] - G(1,1)[n]*G(2,0)[n]);
      pChk = pChk && fabs(detv[k]) > std::numeric_limits<GFTYPE>::epsilon();  // test for zero det
    }
  }
  else {
    for ( n=0; n<detv.size(); n++ ) {
      detv[n] = G(0,0)[n]*(G(1,1)[n]*G(2,2)[n] - G(1,2)[n]*G(2,1)[n])
              - G(0,1)[n]*(G(1,0)[n]*G(2,2)[n] - G(1,2)[n]*G(2,0)[n])
              + G(0,2)[n]*(G(1,0)[n]*G(2,1)[n] - G(1,1)[n]*G(2,0)[n]);
      pChk = pChk && fabs(detv[n]) > std::numeric_limits<GFTYPE>::epsilon();  // test for zero det
    }
  }

} // end of method det


//***********************************************************************************
//***********************************************************************************
// METHOD : Jac_embed
// DESC   : Compute in vectorized way the Jacobian of embedded surface
//          given Gij = dx^i/dxi^j transformation tensor
// ARGS   : 
//          G    : transformation tensor
//          jac  : determinant computed. Must be at least the same length of
//                 each element of G, of size nind, if pind is non-NULLPTR
//          pChk : are all elements positive definite?
//          pind : indirection indices, used if non-NULLPTR
//          nind : number of indirection indices in pind
// RETURNS: none.
//***********************************************************************************
void GElem_base::Jac_embed(GMVFType &G, GTVector<GFTYPE> &jac, GBOOL &pChk, GINT *pind, GINT nind )
{
  GString serr = "GElem_base::Jac_embed: ";

  assert((G.size(1) == 3 && G.size(2) == 3) 
       && "Invalid matrix dimension"); // must have 3 coordinates

  // Compute Jac, and check Jacobian for being positive-definiteness.
  // We compuate the Jacobian as the inner product of
  // (_x_ is vector Cartesian coordinate)
  //      Jacobian is: d_x_/dzeta . _g_
  //
  //                  = G13*gx + G23*gy + G33*gz
  //
  // where,
  //     _g_ = d_x_/dxi X d_x_/deta
  //         = [dy/dxi dz/deta - dz/dxi dy/deta;
  //            dz/dxi dx/deta - dx/dxi dz/deta;
  //            dx/dxi dy/deta - dy/dxi dx/deta]
  //         = [ G21 G32       - G31 G22       ;
  //             G31 G12       - G11 G32       ;
  //             G11 G22       - G21 G12       ]
  // which gives the normal to the embedded surface at all node points
  GSIZET n, k;
  pChk = TRUE;
  GFTYPE xc, yc, zc;
  if ( pind != NULLPTR ) { // 3x3 matrix:
    for ( k=0; k<nind; k++ ) { // loop over desired indices only
      n = pind[k];
      xc = G(1,0)[n]*G(2,1)[n] - G(2,0)[n]*G(1,1)[n];
      yc = G(2,0)[n]*G(0,1)[n] - G(0,0)[n]*G(2,1)[n];
      zc = G(0,0)[n]*G(1,1)[n] - G(1,0)[n]*G(0,1)[n];
      jac[k] = xc*G(0,2)[n] + yc*G(1,2)[n] + zc*G(2,2)[n];
      pChk = pChk && jac[k] > std::numeric_limits<GFTYPE>::epsilon();  // test for zero det
    }
  }
  else {
    for ( n=0; n<jac.size(); n++ ) {
      xc = G(1,0)[n]*G(2,1)[n] - G(2,0)[n]*G(1,1)[n];
      yc = G(2,0)[n]*G(0,1)[n] - G(0,0)[n]*G(2,1)[n];
      zc = G(0,0)[n]*G(1,1)[n] - G(1,0)[n]*G(0,1)[n];
      jac[n] = xc*G(0,2)[n] + yc*G(1,2)[n] + zc*G(2,2)[n];
      pChk = pChk && jac[n] > std::numeric_limits<GFTYPE>::epsilon();  // test for zero det
    }
  }

} // end of method Jac_embed


//***********************************************************************************
//***********************************************************************************
// METHOD : inv
// DESC   : Compute in vectorized way the inverse of specified 
//          (vectorized) matrix
// ARGS   : 
//          G    : matrix to find determinant of. Each element of G is assumed
//                 to have the same length
//          jac  : Jacobian of G (computed prior to entry)
//          iG   : Inverse. Each element must have the same length
// RETURNS: none.
//***********************************************************************************
void GElem_base::inv(GMVFType &G, const GTVector<GFTYPE> &jac, GMVFType &iG)
{
  GString serr = "GElem_base::inv: ";

  assert(((G.size(1) == 2 && G.size(2) == 2) || 
          (G.size(1) == 3 || G.size(2) == 3)) && "Invalid matrix dimension");

  GFTYPE ijac;
  if ( G.size(1) == 2 && G.size(2) == 2 ) {
    for ( GSIZET n=0; n<jac.size(); n++ ) { // 2x2 matrix
      ijac = 1.0/jac[n];
      iG(0,0)[n] =  G(1,1)[n]*ijac;
      iG(0,1)[n] = -G(1,0)[n]*ijac;
      iG(1,0)[n] = -G(0,1)[n]*ijac;
      iG(1,1)[n] =  G(0,0)[n]*ijac;
    }
    return;
  }
  
  for ( GSIZET n=0; n<jac.size(); n++ ) { // 3x3 matrix
    ijac = 1.0/jac[n];
    iG(0,0)[n] =  (G(1,1)[n]*G(2,2)[n]-G(1,2)[n]*G(2,1)[n])*ijac;
    iG(0,1)[n] = -(G(1,0)[n]*G(2,2)[n]-G(2,0)[n]*G(1,2)[n])*ijac;
    iG(0,2)[n] =  (G(1,0)[n]*G(2,1)[n]-G(2,0)[n]*G(1,1)[n])*ijac;
  
    iG(1,0)[n] = -(G(0,1)[n]*G(2,2)[n]-G(2,1)[n]*G(0,2)[n])*ijac;
    iG(1,1)[n] =  (G(0,0)[n]*G(2,1)[n]-G(2,0)[n]*G(0,2)[n])*ijac;
    iG(1,2)[n] = -(G(0,0)[n]*G(2,1)[n]-G(2,0)[n]*G(0,1)[n])*ijac;

    iG(2,0)[n] =  (G(0,1)[n]*G(1,2)[n]-G(1,1)[n]*G(0,2)[n])*ijac;
    iG(2,1)[n] = -(G(0,0)[n]*G(1,2)[n]-G(1,0)[n]*G(0,2)[n])*ijac;
    iG(2,2)[n] =  (G(0,0)[n]*G(1,1)[n]-G(1,0)[n]*G(0,1)[n])*ijac;
  }

} // end of method inv


//***********************************************************************************
//***********************************************************************************
// METHOD : get_indirect
// DESC   : Compute indirection indices
// ARGS   : 
//          b        : basis array
//          vert_ind : indices into data block for vertices
//          edge_ind : indices into data block for edges
//          face_ind : indices into data block for faces
// RETURNS: none.
//***********************************************************************************
void GElem_base::get_indirect(GTVector<GNBasis<GCTYPE,GFTYPE>*> &b, GVVInt &vert_ind,
                             GVVInt &edge_ind, GVVInt &face_ind)
{
  GString serr = "GElem_base::get_indirect: ";

  #if defined(_G_IS1D)
    get_indirect1d(b, vert_ind, edge_ind, face_ind);
  #elif defined(_G_IS2D)
    get_indirect2d(b, vert_ind, edge_ind, face_ind);
  #elif defined(_G_IS2D)
    get_indirect3d(b, vert_ind, edge_ind, face_ind);
  #endif
  
} // end of method get_indirect


//***********************************************************************************
//***********************************************************************************
// METHOD : get_indirect1d
// DESC   : Compute indirection indices for 1d element
// ARGS   : 
//          b        : basis array
//          vert_ind : indices into data block for vertices
//          edge_ind : indices into data block for edges
//          face_ind : indices into data block for faces
// RETURNS: none.
//***********************************************************************************
void GElem_base::get_indirect1d(GTVector<GNBasis<GCTYPE,GFTYPE>*> &b, GVVInt &vert_ind,
                             GVVInt &edge_ind, GVVInt &face_ind)
{
  GString serr = "GElem_base::get_indirect1d: ";

  GTVector<GINT> N(GDIM);
  for ( GSIZET j=0; j<GDIM; j++ ) N[j] = b[j]->getOrder() + 1;

  // Indirection indices for vertices:
  vert_ind.resize(nVertices_);
  for ( GSIZET j=0; j<nVertices_; j++ ) vert_ind[j].resize(1);
  vert_ind[0][0] = 0;
  vert_ind[1][0] = N[0] - 1;

  // Indirection indices for edges:
  edge_ind.resize(nEdges_);
  for ( GSIZET j=0; j<nEdges_; j++ ) edge_ind[j].resize(1);
  edge_ind[0][0] = 0;
  edge_ind[1][0] = N[0] - 1;


  // Indirection indices for faces (same as for edges):
  face_ind.resize(nFaces_);
  for ( GSIZET j=0; j<nFaces_; j++ ) face_ind[j].resize(1);
  face_ind[0][0] = 0;
  face_ind[1][0] = N[0] - 1;
  
} // end of method get_indirect1d


//***********************************************************************************
//***********************************************************************************
// METHOD : get_indirect2d
// DESC   : Compute indirection indices for 2d element
// ARGS   : 
//          b        : basis array
//          vert_ind : indices into data block for vertices
//          edge_ind : indices into data block for edges
//          face_ind : indices into data block for faces
// RETURNS: none.
//***********************************************************************************
void GElem_base::get_indirect2d(GTVector<GNBasis<GCTYPE,GFTYPE>*> &b, GVVInt &vert_ind,
                             GVVInt &edge_ind, GVVInt &face_ind)
{
  GString serr = "GElem_base::get_indirect2d: ";

  GTVector<GINT> N(GDIM);
  for ( GSIZET j=0; j<GDIM; j++ ) N[j] = b[j]->getOrder() + 1;

  // Indirection indices for vertices:
  vert_ind.resize(nVertices_);
  for ( GSIZET j=0; j<nVertices_; j++ ) vert_ind[j].resize(1);
  vert_ind[0][0] = 0;
  vert_ind[1][0] = N[0] - 1;
  vert_ind[2][0] = N[0]*N[1] - 1;
  vert_ind[3][0] = N[0]*(N[1]-1);

  // Indirection indices for edges:
  edge_ind.resize(nEdges_);
  for ( GSIZET j=0; j<nEdges_; j++ ) edge_ind[j].resize(j%2==0?N[0]:N[1]);
  for ( GSIZET j=0; j<N[0]; j++ ) {
    edge_ind[0][j] = j;
    edge_ind[2][j] = N[0]*(N[1]-1) + j;
  }

  for ( GSIZET j=0; j<N[1]; j++ ) {
    edge_ind[1][j] = (j+1)*N[0] - 1;
    edge_ind[3][j] = j*N[0];
  }

  // Indirection indices for faces (same as for edges):
  face_ind.resize(nFaces_);
  for ( GSIZET j=0; j<nFaces_; j++ ) face_ind[j].resize(j%2==0?N[0]:N[1]);
  for ( GSIZET j=0; j<N[0]; j++ ) {
    face_ind[0][j] = j;
    face_ind[2][j] = N[0]*(N[1]-1) + j;
  }

  for ( GSIZET j=0; j<N[1]; j++ ) {
    face_ind[1][j] = (j+1)*N[0] - 1;
    face_ind[3][j] = j*N[0];
  }

  
} // end of method get_indirect2d


//***********************************************************************************
//***********************************************************************************
// METHOD : get_indirect3d
// DESC   : Compute indirection indices for 3d element
// ARGS   : 
//          b        : basis array
//          vert_ind : indices into data block for vertices
//          edge_ind : indices into data block for edges
//          face_ind : indices into data block for faces
// RETURNS: none.
//***********************************************************************************
void GElem_base::get_indirect3d(GTVector<GNBasis<GCTYPE,GFTYPE>*> &b, GVVInt &vert_ind,
                             GVVInt &edge_ind, GVVInt &face_ind)
{
  GString serr = "GElem_base::get_indirect3d: ";

  GTVector<GINT> N(GDIM);
  for ( GSIZET j=0; j<GDIM; j++ ) N[j] = b[j]->getOrder() + 1;

  // Indirection indices for vertices:
  vert_indices_.resize(nVertices_);
  for ( GSIZET j=0; j<nVertices_; j++ ) vert_ind[j].resize(1);
  vert_ind[0][0] = 0;
  vert_ind[1][0] = N[0] - 1;
  vert_ind[2][0] = N[0]*N[1] - 1;
  vert_ind[3][0] = N[0]*(N[1]-1);
  vert_ind[4][0] = N[0]*N[1]*(N[2]-1);
  vert_ind[5][0] = N[0]*N[1]*(N[2]-1) + N[0]-1;
  vert_ind[6][0] = N[0]*N[1]*N[2] - 1;
  vert_ind[7][0] = N[0]*N[1]*N[2] - N_[1];

  // Indirection indices for edges:
  edge_ind.resize(nEdges_);
  for ( GSIZET j=0; j<8; j+=2 ) edge_ind[j].resize(N[0]);
  for ( GSIZET j=0; j<N[0]; j++ ) { // 'x'-horizontal edges
    edge_ind[0][j] = j;
    edge_ind[2][j] = N[0]*(N[1]-1) + j;
    edge_ind[4][j] = N[0]*N[1]*(N[2]-1) + j;
    edge_ind[6][j] = N[0]*N[1]*N[2] - N[0] + j;
  }

  for ( GSIZET j=1; j<8; j+=2 ) edge_ind[j].resize(N[1]);
  for ( GSIZET j=0; j<N[1]; j++ ) { // 'y-'horizontal edges
    edge_ind[1][j] = (j+1)*N[0] - 1;
    edge_ind[3][j] =     j*N[0];
    edge_ind[5][j] = N[0]*N[1]*(N[2]-1) + (j+1)*N[0] - 1;
    edge_ind[7][j] = N[0]*N[1]*(N[2]-1) + j*N[0];
  }

  for ( GSIZET j=8; j<12; j++ ) edge_ind[j].resize(N[2]);
  for ( GSIZET j=0; j<N[2]; j++ ) { // 'z' vertical edges
    edge_ind[8][j] = j    *N[0]*N[1];
    edge_ind[9][j] = j    *N[0]*N[1] + N[0] - 1;
    edge_ind[10][j] = (j+1)*N[0]*N[1] - 1;
    edge_ind[11][j] = (j+1)*N[0]*N[1] - N[0];
  }

  // Indirection indices for faces ...
  face_ind.resize(nFaces_);

  GSIZET i, j, k;

  // ... Face 0 (front):
  face_ind[0].resize(N[0]*N[2]);
  for ( k=0; k<N[2]; k++ )
    for ( i=0; i<N[0]; i++ )
      face_ind[0][i+k*N[0]] = k*N[0]*N[1] + i;

  // ... Face 1 (right side):
  face_ind[1].resize(N[1]*N[2]);
  for ( k=0; k<N[2]; k++ )
    for ( j=0; j<N[1]; j++ )
      face_ind[1][j+k*N[1]] = k*N[0]*N[1] + (j+1)*N[0] - 1;

  // ... Face 2 (back):
  face_ind[2].resize(N[0]*N[2]);
  for ( k=0; k<N[2]; k++ )
    for ( i=0; i<N[0]; i++ )
      face_ind[2][i+k*N[0]] = (k+1)*N[0]*N[1] + i - N[0];

  // ... Face 3 (left side):
  face_ind[3].resize(N[1]*N[2]);
  for ( k=0; k<N[2]; k++ )
    for ( j=0; j<N[1]; j++ )
      face_ind[3][j+k*N[1]] = k*N[0]*N[1] + j*N[0];

  // ... Face 4 (bottom):
  face_ind[4].resize(N[0]*N[1]);
  for ( j=0; j<N[1]; j++ )
    for ( i=0; i<N[0]; i++ )
      face_ind[4][i+j*N[0]] = i + j*N[0];

  // ... Face 5 (top):
  face_ind[5].resize(N[0]*N[1]);
  for ( j=0; j<N[1]; j++ )
    for ( i=0; i<N[0]; i++ )
      face_ind[5][i+j*N[0]] = N[0]*N[1]*(N[2]-1) + i + j*N[1];
  
} // end of method get_indirect3d


//***********************************************************************************
//***********************************************************************************
// METHOD : interp(1)
// DESC   : Do interpolation of specified field that is assumed to reside on
//          this element, to specified reference node locations. Is usually
//          slower.
// ARGS   : 
//          ufrom: basis array to interp to/evauate at. If this is
//                 a 'global' array then the indices gibeg_ and giend_
//                 govern the starting/ending locations for the portion
//                 'owned' by this element, and must be set prior to call. 
//                 The elements of this array (or portion of if) must reside 
//                 on the gbasis_ nodepoints
//          xito : array of same pointers to reference nodes at which to 
//                 interpolate field ufrom
//          uto  : to-vector field at new xito points, in tensor product ordering
// RETURNS: none.
//***********************************************************************************
void GElem_base::interp(GTVector<GFTYPE> &ufrom, GTVector<GTVector<GFTYPE>*> &xito, GTVector<GFTYPE> &uto)
{
  GString serr = "GElem_base::interp(1): ";

  GTVector<GINT> N(GDIM);
  GSIZET i, j, k, n, nnodes;

  for ( j=0,nnodes=1; j<xito.size(); j++ ) nnodes *= xito[j]->size();
  for ( j=0         ; j<xito.size(); j++ ) N[j] = gbasis_[j]->getOrder()+1;

  assert(ufrom.size() < Ntot_ && uto.size() < nnodes && "Inadequate array size(s)");
  gshapefcn_->set_basis(gbasis_);

  Ni_.resize(nnodes);  // shape function resize 
  GTVector<GINT>   I(GDIM);   // tensor produc index

  // NOTE: In general it would be best to compute interpolation
  // using tensor product computations. But this requires tmp arrays
  // that must be passed in depending on whether we are in 1d, 2d or 3d.
  // To reduce the amount of temp data passed, and to prevent additional
  // reallocations, we use the shape function method, and simply sum
  // (_vector_, as opposted to _matrix_) contributions.

  // uto  computed from:
  //    uto_J = Sum_I ufrom^I N_I(xi_J)
  // where N_I(xi_J) is the I-th tensor product shape function 
  // evaluated at xi_J
#if defined(_G_IS1D)
  uto = 0.0;
  for ( i=0, n=0; i<N[0]; i++ ) {
      I[0] = i; 
      gshapefcn_->Ni(I, xito, Ni_); // I-th shape fcn
      uto += Ni_*ufrom[n]; // multiply by old field
  } // i-loop
#endif
#if defined(_G_IS2D)
  uto = 0.0;
  for ( j=0, n=0; j<N[1]; j++ ) {
    for ( i=0; i<N[0]; i++, n++ ) {
      I[0] = i; I[1] = j;
      gshapefcn_->Ni(I, xito, Ni_); // I-th shape fcn
      uto += Ni_*ufrom[n]; // multiply by old field
    } // i-loop
  } // j-loop
#endif
#if defined(_G_IS3D)
  uto = 0.0;
  for ( k=0, n=0; k<N[2]; k++ ) {
    for ( j=0; j<N[1]; j++ ) {
      for ( i=0; i<N[0]; i++, n++ ) {
        I[0] = i; I[1] = j; I[2] = k;
        gshapefcn_->Ni(I, xito, Ni_); // I-th shape fcn
        uto += Ni_*ufrom[n]; // multiply by old field
      } // i-loop
    } // j-loop
  } // k-loop
#endif

  
} // end of method interp(1)


//***********************************************************************************
//***********************************************************************************
// METHOD : interp (2)
// DESC   : Do interpolation of specified field that is assumed to reside on
//          this element, to specified reference node locations. Should usually
//          be faster.
// ARGS   : 
//          ufrom: basis array to interp to/evauate at. If this is
//                 a 'global' array then the indices gibeg_ and giend_
//                 govern the starting/ending locations for the portion
//                 'owned' by this element, and must be set prior to call. 
//                 The elements of this array (or portion of if) must reside 
//                 on the gbasis_ nodepoints
//          xito : array of same pointers to reference nodes at which to 
//                 interpolate field ufrom. If this is a 'global' array, then the
//                 indices gibeg_ and giend_ govern the starting/ending locations
//                 and must be set prior to the call, as for ufrom.
//          uto  : to-vector field at new xito points, in tensor product ordering
//                 If this is a 'global' array, then the indices gibeg_ and giend_ 
//                 govern the starting/ending locations and must be set prior to the call, 
//                 as for ufrom.
//          matv : vector of matrices element of length GDIM. This is temp space, and
//                 the elements will be allocated only when needed
//          matu : vector of matrices element of length GDIM. This is temp space, and
//                 the elements will be allocated only when needed
//          tmp  : temp vector; will be allocated when needed
//          
// RETURNS: none.
//***********************************************************************************
void GElem_base::interp(GTVector<GFTYPE> &ufrom, GTVector<GTVector<GFTYPE>*> &xito, GTVector<GFTYPE> &uto,
GTVector<GTMatrix<GFTYPE>> &matv, GTVector<GTMatrix<GFTYPE>> &matu, GTVector<GFTYPE> &tmp)
{
  GString serr = "GElem_base::interp(2): ";

  GTVector<GINT>   N(GDIM);   // dimensions
  GTVector<GINT>   I(GDIM);   // tensor produc index
  GSIZET j, nnodes;

  for ( j=0,nnodes=1; j<xito.size(); j++ ) nnodes *= xito[j]->size();
  for ( j=0         ; j<xito.size(); j++ ) N[j] = gbasis_[j]->getOrder()+1;

  assert(ufrom.size() < Ntot_ && uto.size() < nnodes && "Inadequate array size(s)");
  assert(matv.size() < GDIM && "Inadequate number of tmp matrices");

#if defined(_G_IS1D)
  matv[0].resize(xito[0]->size(),gbasis_[0]->getOrder()+1);
  gbasis_[0]->evalBasis(*xito[0], matv[0]);
  uto = matv[0]*ufrom;
#endif

#if defined(_G_IS2D)
  for ( j=0; j<GDIM; j++ ) {
    matv[j].resizem(xito[j]->size(),gbasis_[j]->getOrder()+1);
    gbasis_[j]->evalBasis(*xito[j], matv[j]);
  }
  matu[0].resizem(matv[1].size(2),matv[1].size(1));
  matv[1].transpose(matu[0]);
  GMTK::D2_X_D1<GFTYPE>(matv[0], matu[0], ufrom, tmp, uto);
#endif
  
#if defined(_G_IS3D)
  for ( j=0; j<GDIM; j++ ) {
    matv.resizem(xito[j]->size(),gbasis_[j]->getOrder()+1);
    gbasis_[j]->evalBasis(*xito[j], matv[j]);
  }
  matu[0].resizem(matv[1].size(2),matv[1].size(1));
  matu[1].resizem(matv[2].size(2),matv[2].size(1));
  matv[1].transpose(matu[0]); matv[2].transpose(matu[1]);
  GMTK::D3_X_D2_X_D1<GFTYPE>(matv[0], matu[0], matu[1], ufrom, tmp, uto);
#endif

  
} // end of method interp(2)


//***********************************************************************************
//***********************************************************************************
// METHOD : set_faceNormal2d
// DESC   : Compute faceNormal_ components for 2d element
// ARGS   : 
//          rij  : matrix of quantities dx_j/dxi_i
//          
// RETURNS: none.
//***********************************************************************************
void GElem_base::set_faceNormal2d(GTMatrix<GTVector<GFTYPE>> &rij, GTVector<GTVector<GFTYPE>> &faceNormal)
{
  GString serr = "GElem_base::set_faceNormal2d: ";

  // Compute edge-normal vectors at each node edge node point.
  // We compute these as
  //   Normal = d_x_/dxi  X hat(k) = [-dy/dxi ,dx/dxi], sides 0, 2
  //   Normal = d_x_/deta X hat(k) = [-dy/deta,dx/deta] sides 1, 3
  // when not embedded, and
  //   Normal = d_x_/dxi  X _x_, sides 0, 2
  //   Normal = d_x_/deta X _x_, sides 1, 3
  // if embedded (and normalize; _x_ is vector Cartesian coord)

  GSIZET nxy = elemtype_ == GE_2DEMBEDDED ? GDIM+1 : GDIM;
  GSIZET istart;

  GFTYPE *fn[2*GDIM];   // faceNormal data pointers
  for ( GSIZET j=0; j<nxy; j++ ) fn[j] = faceNormal[j].data();

  if ( elemtype_ == GE_2DEMBEDDED ) { // embedded surfaces

    // These must be computed in the given order:
    istart = 0;
    GMTK::cross_prod<GFTYPE>(rij         (0,0), rij         (1,0), rij             (2,0),
                             xNodes_       [0], xNodes_       [1], xNodes_           [2], 
                             edge_indices_ [0].data()            , edge_indices_[0].size(), 
                             fn[0]+istart     , fn[1]+istart     , fn[2]+istart           );
    istart += edge_indices_[0].size();
    GMTK::cross_prod<GFTYPE>(rij         (0,1), rij         (1,1), rij             (2,1),
                             xNodes_       [0], xNodes_       [1], xNodes_           [2], 
                             edge_indices_ [1].data()            , edge_indices_[1].size(), 
                             fn[0]+istart     , fn[1]+istart     , fn[2]+istart           );

    istart += edge_indices_[1].size();
    GMTK::cross_prod<GFTYPE>(rij         (0,0), rij         (1,0), rij             (2,0),
                             xNodes_       [0], xNodes_       [1], xNodes_           [2], 
                             edge_indices_ [2].data()            , edge_indices_[2].size(), 
                             fn[0]+istart     , fn[1]+istart     , fn[2]+istart           );

    istart += edge_indices_[2].size();
    GMTK::cross_prod<GFTYPE>(rij         (0,1), rij         (1,1), rij             (2,1),
                             xNodes_       [0], xNodes_       [1], xNodes_           [2], 
                             edge_indices_ [3].data()            , edge_indices_[3].size(), 
                             fn[0]+istart     , fn[1]+istart     , fn[2]+istart           );
  }
  else { // non-embedded surfaces

    // These must be computed in the given order:
    istart  = 0;
    GMTK::cross_prod_k<GFTYPE>(rij         (0,0), rij         (1,0), 
                               edge_indices_ [0].data() , edge_indices_[0].size(),  1,
                               fn[0]+istart     , fn[1]+istart     ); 

    istart += edge_indices_[0].size();
    GMTK::cross_prod_k<GFTYPE>(rij         (0,1), rij         (1,1), 
                               edge_indices_ [1].data() , edge_indices_[1].size(),  1,
                               fn[0]+istart     , fn[1]+istart     ); 

    istart += edge_indices_[1].size();
    GMTK::cross_prod_k<GFTYPE>(rij         (0,0), rij         (1,0), 
                               edge_indices_ [2].data() , edge_indices_[2].size(), -1,
                               fn[0]+istart     , fn[1]+istart     ); 


    istart += edge_indices_[2].size();
    GMTK::cross_prod_k<GFTYPE>(rij         (0,1), rij         (1,1), 
                               edge_indices_ [3].data() , edge_indices_[3].size(), -1,
                               fn[0]+istart     , fn[1]+istart     ); 
  }
  
} // end of method set_faceNormal2d


//***********************************************************************************
//***********************************************************************************
// METHOD : set_faceNormal3d
// DESC   : Compute faceNormal  components for 3d element
// ARGS   : 
//          rij  : matrix of quantities dx_j/dxi_i
//          
// RETURNS: none.
//***********************************************************************************
void GElem_base::set_faceNormal3d(GTMatrix<GTVector<GFTYPE>> &rij, GTVector<GTVector<GFTYPE>> &faceNormal)
{
  GString serr = "GElem_base::set_faceNormal3d: ";

  GSIZET istart;

  GFTYPE *fn[6]; // pointers to faceNormal data
  for ( GSIZET k=0; k<GDIM; k++ ) fn[k] = faceNormal[k].data();

  // For each face, take cross prod of d_x_/dxi_1 X d_x_/dxi_2
  // where xi1, and xi2 are the reference coords of the face:
  istart = 0; // where to start in faceNormal; order matters!
  for ( GSIZET k=0; k<4; k++ ) { // vertical faces
    GMTK::cross_prod<GFTYPE>(rij           (0,k%2), rij         (1,k%2), rij             (2,k%2),
                             rij             (0,2), rij           (1,2), rij               (2,2),
                             face_indices_ [k].data()                  , face_indices_[k].size(), 
                             fn[0]+istart         , fn[1]+istart       , fn[2]+istart          );
    istart += face_indices_[k].size();
  }
  
  for ( GSIZET k=4; k<6; k++ ) { // bottom, top faces
    GMTK::cross_prod<GFTYPE>(rij             (0,0), rij           (1,0), rij               (2,0),
                             rij             (0,1), rij           (1,1), rij               (2,1),
                             face_indices_ [k].data()                  , face_indices_[k].size(), 
                             fn[0]+istart         , fn[1]+istart       , fn[2]+istart          );
    istart += face_indices_[k].size();
  }
  
} // end of method set_faceNormal3d


//***********************************************************************************
//***********************************************************************************
// METHOD : init
// DESC   : Initialize element. set_basis, set_elemtype should already have 
//          been called if necessary prior to entry
// ARGS   : 
//          xNodes: Cartesian coordinates for all node points (each
//                  node point is in tensor product form).
// RETURNS: none.
//***********************************************************************************
void GElem_base::init(GVVFType &xNodes)
{
  GString serr = "GElem_base::init: ";

  assert( bbasis_ && elemtype_ != GE_MAX && "Basis or element type not set");
  
  // Set xNodes_ data member. Note that the input array sizes may be
  // larger than for xNodes_:
  for ( GSIZET j=0; j<xNodes.size(); j++ ) {
    for ( GSIZET i=0; i<xNodes_[j].size(); i++ ) xNodes_[j][i] = xNodes[j][i];
  }
  
  build_elem();
  bInitialized_ = TRUE;
  
} // end of method init


//**********************************************************************************
//**********************************************************************************
// METHOD : Assignment method
// DESC   : 
// ARGS   : GLLBasis &
// RETURNS: none
//**********************************************************************************
void GElem_base::operator=(const GElem_base &e)
{
  if ( &e == this ) return;
 
 // copy data:
  Ntot_         = e.Ntot_;
  N_            = e.N_;
  igbeg_        = e.igbeg_;
  igend_        = e.igend_;
  bInitialized_ = e.bInitialized_;
  elemtype_     = e.elemtype_;
  elemid_       = e.elemid_;
  rootid_       = e.rootid_;
  nVertices_    = e.nVertices_;
  nEdges_       = e.nEdges_;
  nFaces_       = e.nFaces_;
  gbasis_       = e.gbasis_;
  xNodes_       = e.xNodes_;
  xiNodes_      = e.xiNodes_;
  xVertices_    = e.xVertices_;
  volume_       = e.volume_;
  gshapefcn_    = e.gshapefcn_;
  edgeCentroid_ = e.edgeCentroid_;
  faceCentroid_ = e.faceCentroid_;
  elemCentroid_ = e.elemCentroid_;
#if 0
  dXidX_        = e.dXidX_;
  Jac_          = e.Jac_;
  faceJac_      = e.faceJac_;
  faceNormal_   = e.faceNormal_;
#endif
  vert_indices_ = e.vert_indices_;
  edge_indices_ = e.edge_indices_;
  face_indices_ = e.face_indices_;
  bdy_indices_  = e.bdy_indices_;
  bdy_types_    = e.bdy_types_;
  mask_         = e.mask_;

} // end, method operator=
