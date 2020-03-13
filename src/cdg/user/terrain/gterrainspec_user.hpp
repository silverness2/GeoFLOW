//==================================================================================
// Module       : gterrainspec.hpp
// Date         : 7/11/19 (DLR)
// Description  : Terrain specification function implementations provided
//                by user
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GTERRAINSPEC_USER_HPP)
#define _GTERRAINSPEC_USER_HPP

#include "tbox/property_tree.hpp"
#include "gtypes.h"
#include "gtvector.hpp"
#include "ggrid.hpp"


using namespace geoflow;
using namespace geoflow::tbox;

typedef GFTYPE                      Time;
typedef GTVector<GTVector<GFTYPE>*> State;


namespace gterrainspec
{

GBOOL impl_gauss_range     (const PropertyTree& stree, GGrid &grid,  State &utmp, State &xb);
GBOOL impl_poly_range      (const PropertyTree& stree, GGrid &grid,  State &utmp, State &xb);
GBOOL impl_schar_range     (const PropertyTree& stree, GGrid &grid,  State &utmp, State &xb);
};

#endif
