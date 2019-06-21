//==================================================================================
// Module       : ggrid_box
// Date         : 11/11/18 (DLR)
// Description  : Object defining a (global) 2d or 3d box grid.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#if !defined(_GGRID_BOX_HPP)
#define _GGRID_BOX_HPP

#include "gtypes.h"
#include <functional>
#include "gtvector.hpp"
#include "gtmatrix.hpp"
#include "gnbasis.hpp"
#include "gelem_base.hpp"
#include "gdd_base.hpp"
#include "ggrid.hpp"
#include "gshapefcn_linear.hpp"
#include "polygon.h"
#include "gtpoint.hpp"

typedef GTMatrix<GFTYPE> GFTMatrix;

class GGridBox 
{

public:
                            GGridBox(GTPoint<GFTYPE> &P0, GTPoint<GFTYPE> &P1, GTVector<GINT> &ne, GTVector<GNBasis<GCTYPE,GFTYPE>*> &b, GINT nprocs); // 2d, 3d constructor
                           ~GGridBox();

        void                do_grid(GGrid &grid, GINT irank);             // compute grid for irank
        void                set_partitioner(GDD_base *d);                 // set and use GDD object
        void                set_bdy_callback(
                            std::function<void(GGrid &)> &callback);     // set bdy-set callback
        void                set_basis(GTVector<GNBasis<GCTYPE,GFTYPE>*> &b);             // set element basis

        void                print(const GString &filename);              // print grid to file

friend  std::ostream&       operator<<(std::ostream&, GGridBox &);       // Output stream operator
 

//protected:
         void               init2d();                                       // initialize base icosahedron for 2d grid
         void               init3d();                                       // initialize for 3d grid
         void               reorderverts2d(GTVector<GTPoint<GFTYPE>> &, GTVector<GSIZET>&,
                                           GTVector<GTPoint<GFTYPE>> &);    // make verts consis with shapefcns


         void               do_grid2d(GGrid &grid, GINT rank);            // do 2d grid
         void               do_grid3d(GGrid &grid, GINT rank);            // do 3d grid
       

private:

GINT                    ndim_;          // grid dimensionality (2 or 3)
GINT                    nprocs_;        // no. MPI tasks
GDD_base               *gdd_;           // domain decomposition/partitioning object
GShapeFcn_linear       *lshapefcn_;     // linear shape func to compute 2d coords
GTPoint<GFTYPE>         P0_;            // P0 = starting point of box origin
GTPoint<GFTYPE>         P1_;            // P1 = diagonally-opposing box point
GTVector<GTPoint<GFTYPE>>
                        ftcentroids_;   // centroids of finest elements
GTVector<GNBasis<GCTYPE,GFTYPE>*> 
                        gbasis_;        // directional bases
GTVector<GQuad<GFTYPE>> qmesh_;         // list of vertices for each 2d (quad) element
GTVector<GHex<GFTYPE>>  hmesh_;         // list of vertices for each 3d (hex) element
GTVector<GINT>          ne_;            // # elems in each coord direction in 3d
GTVector<GFTYPE>        Lbox_;          // length of box edges (x, y, [and z])
std::function<void(GGrid&)>
                       *bdycallback_;   // callback object+method to set bdy conditions

};

#endif