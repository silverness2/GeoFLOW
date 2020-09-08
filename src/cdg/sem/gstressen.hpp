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
        using Ftype      = typename Interface::Value;
        using Derivative = typename Interface::Derivative;
        using Time       = typename Interface::Time;
        using CompDesc   = typename Interface::CompDesc;
        using Jacobian   = typename Interface::Jacobian;
        using Size       = typename Interface::Size;

        static_assert(std::is_same<State,GTVector<GTVector<Ftype>*>>::value,
               "State is of incorrect type");
        static_assert(std::is_same<StateComp,GTVector<Ftype>>::value,
               "StateComp is of incorrect type");
        static_assert(std::is_same<Derivative,GTVector<GTVector<Ftype>*>>::value,
               "Derivative is of incorrect type");
        static_assert(std::is_same<Grid,GGrid>::value,
               "Grid is of incorrect type");

public:

                          GStressEnOp() = delete;
                          GStressEnOp(Grid &grid);
                          GStressEnOp(const GStressEnOp &);
                         ~GStressEnOp();

        void              apply(State &u, GINT idir, State  &utmp, 
                                StateComp &si);                              // stress op evaluation in idir
        void              apply(State &u, State  &utmp,  
                                StateComp &e);                               // stress-energy op evaluation
        void              set_mu(StateComp &mu);                             // set viscosity
        void              set_kappa(StateComp &kappa);                       // set energy kappa


private:
        GBOOL                        bown_mu_;    // flag telling instance if it owns mu_
        GBOOL                        bown_kappa_; // flag telling instance if it owns kappa_
        Mass                         *massop_;    // mass matrix, required
        Grid                         *grid_;      // grid set on construction
        StateComp                    *mu_;        // viscosity
        StateComp                    *kappa_;     // viscosity


};


#include "gstressen.ipp"


#endif