//==================================================================================
// Module       : gcg.hpp
// Date         : 3/7/20 (DLR)
// Description  : Encapsulates the methods and data associated with
//                a conjugate Gradient (CG) solver
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================

#if !defined(_GCG_HPP)
#define _GCG_HPP

#include <sstream>
#include "gtvector.hpp"
#include "gtmatrix.hpp"
#include "glinop.hpp"
#include "gcomm.hpp"
#include "ggfx.hpp"


template<typename TypePack> 
class GCG : public LinSolverBase<TypePack>
{
public:
                      enum  GCGERR         {GCGERR_NONE=0, GCGERR_NOCONVERGE, GCGERR_SOLVE, GCGERR_PRECOND,  GCGERR_BADITERNUM};
                      using Types          = TypePack;
                      using Operator       = typename Types::Operator;
                      using Preconditioner = typename Types::Preconditioner;
                      using State          = typename Types::State;
                      using StateComp      = typename Types::StateComp;
                      using Grid           = typename Types::Grid;
                      using Value          = typename Types::Value;
                      using ConnectivityOp = typename Types::ConnectivityOp;

                      static_assert(std::is_same<Operator,GLinOp>::value,
                                    "Operator is of incorrect type");
                      static_assert(std::is_same<State,GTVector<GTVector<GFTYPE>*>>::value,
                                    "State is of incorrect type");
                      static_assert(std::is_same<StateComp,GTVector<GFTYPE>>::value,
                                    "StateComp is of incorrect type");
                      static_assert(std::is_same<Value,GFTYPE>::value,
                                    "Value is of incorrect type");
                      static_assert(std::is_same<ConnectivityOp,GGFX<Value>>::value,
                                    "ConnectivityOp is of incorrect type");



                      GCG() = delete;
                      GCG(Traits& traits, Grid& grid, ConnectivityOp& ggfx, 
                          State& tmppack);
                     ~GCG();
                      GCG(const GCG &a);
                      GCG  &operator=(const GCG &);
                   

                      void         solve_impl(Operator& A, const StateComp& b, StateComp& x);
                      void         solve_impl(Operator& A, const StateComp& b, 
                                              const StateComp& xb, StateComp& x);
                      void         solve_impl(const StateComp& b, StateComp& x) {};
                      StateComp&   get_residuals() { return residuals_; }  
                      GINT         get_iteration_count() { return iter_+1; }  


private:
// Private methods:
                      void     init();
                      GFTYPE   compute_norm(const StateComp&, State&);
// Private data:
     GC_COMM           comm_;        // communicator
     GINT              irank_;       // rank
     GINT              iter_;        // iteration number
     GINT              nprocs_;      // no. tasks
     StateComp         residuals_;   // list of resituals for each iteration
     LinSolverBase<TypePack>
                      *precond_;     // preconditioner


};

#include "gcg.ipp"

#endif

