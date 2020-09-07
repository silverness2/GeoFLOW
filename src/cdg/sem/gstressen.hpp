//==================================================================================
// Module       : gstress.hpp
// Date         : 09/05/20 (DLR)
// Description  : Represents the SEM discretization of the full viscous
//                stress-energy operator. The viscous stress in the 
//                momentum eqution is
//                    2 [ mu s_{ij}],j,
//                where
//                    s_{ij} = (Del_i u_j + Del_j u_i)/2
//                and the viscous stress-energy for the energy equation is
//                    2 [mu u_i s_{ij} ],j
//                where u_i is the velocity, and mu, the viscosity. Repeated
//                indices are summed here.  For the energy, this is a nonlinear 
//                operator, so should not derive from GLinOp. Operator requires 
//                that grid consist of elements of only one type.
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none
//==================================================================================

#if !defined(_GSTRESSENERGYOP_HPP)
#define _GSTRESSENERGYOP_HPP
#include "gtvector.hpp"
#include "gnbasis.hpp"
#include "ggrid.hpp"
#include "gmass.hpp"
#include "gtmatrix.hpp"
#include "gmtk.hpp"
#include "pdeint/equation_base.hpp"


template<typename TypePack>
class GStressEnOp
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

                          GstressEnOp() = delete;
                          GstressEnOp(Grid &grid);
                          GstressEnOp(const GstressEnOp &);
                         ~GstressEnOp();

        void              apply(State &u, GINT idir, State  &utmp, 
                                StateComp &si);                              // stress op evaluation in idir
        void              apply(State &u, State  &utmp,  
                                StateComp &e);                               // stress-energy op evaluation
        void              set_mu(StateComp &mu);                             // set viscosity


private:
        Mass                         *massop_; // mass matrix, required
        Grid                         *grid_;   // grid set on construction
        StateComp                    *mu_;     // viscosity


};


#include "gstressen.ipp"


#endif
