//==================================================================================
// Module       : gterrainspec.hpp
// Date         : 7/11/19 (DLR)
// Description  : Terrain specification function implementations provided
//                by user for box grids (2d and 3d).
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GTERRAINSPECBOX_USER_HPP)
#define _GTERRAINSPECBOX_USER_HPP

#include "tbox/property_tree.hpp"
#include "gtypes.h"
#include "gtvector.hpp"
#include "ggrid.hpp"


using namespace geoflow;
using namespace geoflow::tbox;

typedef GFTYPE                      Time;
typedef GTVector<GTVector<GFTYPE>*> State;


namespace gterrain_specbox
{

GBOOL impl_gauss_range     (const PropertyTree& ptree, GGrid &grid,  State &utmp, State &xb);
GBOOL impl_poly_range      (const PropertyTree& ptree, GGrid &grid,  State &utmp, State &xb);
GBOOL impl_schar_range     (const PropertyTree& ptree, GGrid &grid,  State &utmp, State &xb);
};

#endif
