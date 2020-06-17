//==================================================================================
// Module       : gpdv.hpp
// Date         : 11/11/18 (DLR)
// Description  : Represents the SEM discretization of the 'pdV' operator:
//                p Div u. This is a nonlinear operator, so should not derive 
//                from GLinOp. This operator requires that grid consist of
//                elements of only one type.
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none
//==================================================================================

#if !defined(_GLAPLACIANOP_HPP)
#define _GLAPLACIANOP_HPP
#include "gtvector.hpp"
#include "gnbasis.hpp"
#include "ggrid.hpp"
#include "gmass.hpp"
#include "gtmatrix.hpp"
#include "gmtk.hpp"


template<typename TypePack>
classGpdv 
{
public:
        using Interface  = EquationBase<TypePack>;
        using Base       = Interface;
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
        static_assert(std::is_same<StateComp,GTVector<Value>::value,
               "StateComp is of incorrect type");
        static_assert(std::is_same<Derivative,GTVector<GTVector<Value>*>>::value,
               "Derivative is of incorrect type");
        static_assert(std::is_same<Grid,GGrid>::value,
               "Grid is of incorrect type");

public:

                          GpdV(Grid &grid, Mass &massop);
                          GpdV(const GpdV &);
                         ~GpdV();

        void              apply(StateComp &p, State &u, State  &utmp, 
                                StateComp &po);                              // Operator-field evaluation
        void              init();                                            // must call after all 'sets'

private:
        void              def_init();
        void              reg_init();
        void              def_prod(StateComp &p State &u, State *utmp,
                                   StateComp &po);
        void              reg_prod(StateComp &p, State &u,  State &utmp, 
                                   StateComp &po);

        GBOOL                        bInitialized_;
        StateComp                    etmp1_;  // elem-based (non-global) tmp vector
        State                        G_;      // metric components
        Mass                         *massop_; // mass matrix, required
        Grid                         *grid_;   // grid set on construction


};


#include "gpdV.ipp"


#endif
