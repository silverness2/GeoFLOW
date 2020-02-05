//==================================================================================
// Module       : gposixio_observer.hpp
// Date         : 3/18/19 (DLR)
// Description  : Observer object for carrying out binary output of state.
// Copyright    : Copyright 2019. Colorado State University. All rights reserved.
// Derived From : ObserverBase.
//==================================================================================
#if !defined(_GPOSIXIO_OBSERVER_HPP)
#define _GPOSIXIO_OBSERVER_HPP

#include "gtvector.hpp"
#include "ggrid.hpp"
#include "pdeint/equation_base.hpp"
#include "pdeint/observer_base.hpp"
#include "tbox/property_tree.hpp"
#include "tbox/mpixx.hpp"
#include "gio.h"

using namespace geoflow::pdeint;
using namespace std;


template<typename EquationType>
class GIOObserver : public ObserverBase<EquationType>
{

public:
        using Equation    = EquationType;
        using EqnBase     = EquationBase<EquationType>;
        using EqnBasePtr  = std::shared_ptr<EqnBase>;
        using IOBaseType  = IOBase<EquationType>;
        using IOBaseTraits= IOBase<EquationType>::Traits;
        using IOBasePtr   = std::shared_ptr<IOBaseType>;
        using State       = typename Equation::State;
        using StateInfo   = typename Equation::StateInfo;
        using Grid        = typename Equation::Grid;
        using Value       = typename Equation::Value;
        using Derivative  = typename Equation::Derivative;
        using Time        = typename Equation::Time;
        using Jacobian    = typename Equation::Jacobian;
        using Size        = typename Equation::Size;

//      using ObserverBase<EquationType>::ObsType;
//      using OBS_CYCLE = typename ObserverBase<EquationType>::ObsType::OBS_CYCLE;
//      using OBS_TIME  = typename ObserverBase<EquationType>::OBS_TIME;

        static_assert(std::is_same<State,GTVector<GTVector<GFTYPE>*>>::value,
               "State is of incorrect type");
        static_assert(std::is_same<Derivative,GTVector<GTVector<GFTYPE>*>>::value,
               "Derivative is of incorrect type");
        static_assert(std::is_same<Grid,GGrid>::value,
               "Grid is of incorrect type");

                           GIOObserver() = delete;
                           GIOObserver(const EqnBasePtr &equation, const IOBasePtr &io_base_ptr,
                                            Grid &grid, typename ObserverBase<EquationType>::Traits &traits);
                          ~GIOObserver() = default;
                           GIOObserver(const GIOObserver &a) = default;
                           GIOObserver &operator=(const GIOObserver &bu) = default;

        void               observe_impl(const Time &t, const State &u, const State &uf);

        void               init(StateInfo &info);
private:
// Private methods:
        void               print_derived(const Time &t, const State &u, GIOTraits &traits, const GC_COMM &comm);
// Private data:
        GBOOL                bprgrid_       ;// print grid flag
        GBOOL                bInit_         ;// is initialized?
        GSIZET               cycle_last_    ;// most recent output cycle
        GSIZET               cycle_         ;// continuously-running cycle
        GSIZET               ocycle_        ;// output cycle number
        GFTYPE               time_last_     ;// most recent output time
        std::vector<GINT>    def_stateindex_;// list of state indices to print
        std::vector<GString> def_statenames_;// list of names of state files
        std::vector<GString> def_gridnames_ ;// list of names of grid comp files
        State                up_            ;// print-state
        State                gp_            ;// print-grid 
        IOBasePtr           *io_ptr_        ;// ptr to IO object
//      IOBaseTraits        *iotraits_      ;// IOBase traits pointer
        StateInfo            stateinfo_     ;// info structure
        StateInfo            grstateinfo_   ;// info structure

};

#include "gposixio_observer.ipp"

#endif

