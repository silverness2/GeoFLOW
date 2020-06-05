//==================================================================================
// Module       : ginitfv.hpp
// Date         : 7/16/19 (DLR)
// Description  : Velocity force inititial condition implementations 
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================
#if !defined(_GINITFV_HPP)
#define _GINITFV_HPP

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


namespace ginitfv
{

GBOOL impl_rand      (const PropertyTree &ptree, GString &sconfig, GGrid &grid, StateInfo &stinfo, Time &time, State &utmp, State &ub, State &uf);
GBOOL impl_abc_box   (const PropertyTree &ptree, GString &sconfig, GGrid &grid, StateInfo &stinfo, Time &time, State &utmp, State &ub, State &uf);
GBOOL impl_abc_icos  (const PropertyTree &ptree, GString &sconfig, GGrid &grid, StateInfo &stinfo, Time &time, State &utmp, State &ub, State &uf);

};



#endif
