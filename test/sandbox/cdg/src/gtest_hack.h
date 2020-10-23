//==================================================================================
// Module       : gtest_hack.cpp
// Date         : 10/21/20 (DLR)
// Description  : GeoFLOW mini-app for tensor-producs performance 
//                improvement
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : 
//==================================================================================

#include "gexec.h"
#include "gtypes.h"
#include <cstdio>
#include <unistd.h>
#include <iostream>
#if defined(_G_USE_GPTL)
  #include "gptl.h"
#endif
#include <memory>
#include <cstdlib>
#include <cassert>
#include <random>
#include "gcomm.hpp"
#include "ggfx.hpp"
#include "gllbasis.hpp"
#include "gmass.hpp"
#include "gmorton_keygen.hpp"
#include "gmtk.hpp"
#include "ggrid_factory.hpp"
#include "pdeint/observer_base.hpp"
#include "pdeint/observer_factory.hpp"
#include "tbox/property_tree.hpp"
#include "tbox/mpixx.hpp"
#include "tbox/global_manager.hpp"
#include "tbox/input_manager.hpp"
#include "gtools.h"


#if 0
typedef  GTVector<GTVector<GFTYPE>*> State;
typedef  GTVector<GFTYPE>            StateComp;
#endif

using namespace geoflow::tbox;
using namespace std;



template< // Complete typepack
typename StateType     = GTVector<GTVector<GFTYPE>*>,
typename StateCompType = GTVector<GFTYPE>,
typename StateInfoType = GStateInfo,
typename GridType      = GGrid,
typename MassOpType    = GMass,
typename ValueType     = GFTYPE,
typename DerivType     = StateType,
typename TimeType      = ValueType,
typename CompType      = GTVector<GStateCompType>,
typename JacoType      = StateType,
typename SizeType      = GSIZET
>
struct TypePack {
        using State      = StateType;
        using StateComp  = StateCompType;
        using StateInfo  = StateInfoType;
        using Grid       = GridType;
        using Mass       = MassOpType;
        using Value      = ValueType;
        using Derivative = DerivType;
        using Time       = TimeType;
        using CompDesc   = CompType;
        using Jacobian   = JacoType;
        using Size       = SizeType;
};
using MyTypes       = TypePack<>;           // Define grid types used
using Grid          = GGrid;
using IOBaseType    = IOBase<MyTypes>;          // IO Base type
using IOBasePtr     = std::shared_ptr<IOBaseType>;// IO Base ptr
using ObsTraitsType = ObserverBase<MyTypes>::Traits;

Grid      *grid_ = NULLPTR;
IOBasePtr  pIO_  = NULLPTR; // ptr to IOBase operator

