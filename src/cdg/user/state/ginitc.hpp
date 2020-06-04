//==================================================================================
// Module       : ginitc.hpp
// Date         : 7/16/19 (DLR)
// Description  : Prescribed (vec field) inititial condition implementations 
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GINITC_HPP)
#define _GINITC_HPP

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



namespace ginitc
{

GBOOL impl_rand      (const PropertyTree &ptree, GString &sconfig, GGrid &grid, StateInfo &stinfo, Time &time, State &utmp, State &ub, State &u);

};



#endif
