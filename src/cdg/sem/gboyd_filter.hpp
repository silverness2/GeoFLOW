//==================================================================================
// Module       : gboyd_filter.hpp
// Date         : 9/14/20 (DLR)
// Description  : Computes the Boyd filter to diminish aliasing errors.
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : FilterBase
//==================================================================================

#if !defined(_GBOYDFILTER_HPP)
#define _GBOYDFILTER_HPP
#include "gtvector.hpp"
#include "gnbasis.hpp"
#include "ggrid.hpp"
#include "gmass.hpp"
#include "gtmatrix.hpp"
#include "gmtk.hpp"
#include "pdeint/filter_base.hpp"


template<typename TypePack>
class GBoydFilter : public FilterBase
{
public:
        using Interface  = EquationBase<TypePack>;
        using State      = typename Interface::State;
        using StateComp  = typename Interface::StateComp;
        using Grid       = typename Interface::Grid;
        using Mass       = typename Interface::Mass;
        using Value      = typename Interface::Value;
        using Derivative = typename Interface::Derivative;
        using Time       = typename Interface::Time;
        using CompDesc   = typename Interface::CompDesc;
        using Jacobian   = typename Interface::Jacobian;
        using Size       = typename Interface::Size;

        static_assert(std::is_same<State,GTVector<GTVector<Value>*>>::value,
               "State is of incorrect type");
        static_assert(std::is_same<StateComp,GTVector<Value>>::value,
               "StateComp is of incorrect type");
        static_assert(std::is_same<Derivative,GTVector<GTVector<Value>*>>::value,
               "Derivative is of incorrect type");
        static_assert(std::is_same<Grid,GGrid>::value,
               "Grid is of incorrect type");

public:

                          GBoydFilter() = delete;
                          GBoydFilter(Grid &grid);
                          GBoydFilter(const GBoydFilter &);
                         ~GBoydFilter();

        void              apply(Time t, StateComp &u, State  &utmp, 
                                StateComp &po);
        void              apply(Time t, StateComp &u, State  &utmp); 

private:
        void              init();

        GBOOL                         bInit_;    // is filter initialized?
        GINT                          ifilter_;  // filter mode
        Ftype                         mufilter_; // truncation amount 
        GTMatrix<Ftype>               Lambda_;   // mode-weighting matrix
        Grid                         *grid_;     // grid set on construction


};


#include "gboyd_filter.ipp"


#endif
