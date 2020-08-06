//==================================================================================
// Module       : ginits.hpp
// Date         : 7/16/19 (DLR)
// Description  : Scalar inititial condition implementations 
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GINITS_HPP)
#define _GINITS_HPP

#include "tbox/property_tree.hpp"
#include "gtypes.h"
#include "gstateinfo.hpp"
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "gutils.hpp"


using namespace geoflow;
using namespace geoflow::tbox;

typedef GFTYPE                      Time;
typedef GTVector<GTVector<GFTYPE>*> State;
typedef GStateInfo                  StateInfo;



namespace ginits
{

GBOOL impl_rand      (const PropertyTree &ptree, GString &sconfig, GGrid &grid, StateInfo &stinfo, Time &time, State &utmp, State &ub, State &u);

};



#endif
