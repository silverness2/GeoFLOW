//==================================================================================
// Module       : gmconv.hpp
// Date         : 6/11/20 (DLR)
// Description  : Object defining a moist convection solver:
//
//                PDEs:
//                     TBD
//
//                This solver can be built in 2D or 3D for box grids,
//                but is valid only for 3D spherical grids.
//
//                The State vector consists of the following:
//                  vx/sx   : x-velocity or momentum density 
//                  vy/sy   : y-velocity or momentum density 
//                  vz/sz   : z-velocity or momentum density 
//                  e       : internal or total energy density
//                  rho_tot : total density or dry density, if no moisture
//                  qvapor  : water vapor mass fraction
//                  q_liq_0 : liquid substance 0 mass fraction |
//                  q_liq_1 : liquid substance 1 mass fraction |  'liquid' mass sector
//                  q_liq_2 : liquid substance 2 mass fraction |
//                   ...
//                  q_ice_0 : 'ice' substance 0 mass fraction  |
//                  q_ice_1 : 'ice' substance 1 mass fraction  |  'ice' mass sector
//                  q_ice_2 : 'ice' substance 2 mass fraction  |
//                   ...
//                  rho_s   : base state density               |
//                  p_s     : base state pressure              | base state sector
//                  w_liq_0 : liquid substance 2 term velocity |
//                  w_liq_1 : liquid substance 2 term velocity | 'liquid' term vel. sector
//                  w_liq_2 : liquid substance 2 term velocity |

//                  w_ice_0 : 'ice' substance 2 term velocity  |
//                  w_ice_1 : liquid substance 2 term velocity | 'ice' term vel. sector
//                  w_ice_2 : 'ice' substance 2 term velocity  |
//                   ...
//
//                The terminal velocities in this state are prescribed,
//                if used; all others are evolved. If hydrometeor fallout
//                is specified, then terminal velocities for all hydrometeors
//                must be provided
// 
// 
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : EquationBase.
//==================================================================================
#if !defined(_GMCONV_HPP)
#define _GMCONV_HPP

#include "gtypes.h"
#include "gphysics.h"
#include <functional>
#include <cstdlib>
#include <memory>
#include <cmath>
#include "gtvector.hpp"
#include "gdd_base.hpp"
#include "ggrid.hpp"
#include "gab.hpp"
#include "gext.hpp"
#include "gbdf.hpp"
#include "gpdv.hpp"
#include "gmass.hpp"
#include "gadvect.hpp"
#include "ghelmholtz.hpp"
//#include "gflux.hpp"
#include "gexrk_stepper.hpp"
#include "gbutcherrk.hpp"
#include "ggrid_box.hpp"
#include "ggrid_icos.hpp"
#include "gmtk.hpp"
#include "ggfx.hpp"
#include "pdeint/equation_base.hpp"

using namespace geoflow::pdeint;
using namespace std;

template<typename TypePack>
class GMConv : public EquationBase<TypePack>
{
public:
        using Interface  = EquationBase<TypePack>;
        using Base       = Interface;
        using State      = typename Interface::State;
        using StateComp  = typename Interface::StateComp;
        using Grid       = typename Interface::Grid;
        using Ftype      = typename Interface::Value;
        using Derivative = typename Interface::Derivative;
        using Time       = typename Interface::Time;
        using CompDesc   = typename Interface::CompDesc;
        using Jacobian   = typename Interface::Jacobian;
        using Size       = typename Interface::Size;

        static_assert(std::is_same<State,GTVector<GTVector<GFTYPE>*>>::value,
               "State is of incorrect type");
        static_assert(std::is_same<StateComp,GTVector<GFTYPE>>::value,
               "StatCompe is of incorrect type");
        static_assert(std::is_same<Derivative,GTVector<GTVector<GFTYPE>*>>::value,
               "Derivative is of incorrect type");
        static_assert(std::is_same<Grid,GGrid>::value,
               "Grid is of incorrect type");

        // MConv solver traits:
        struct Traits {
          GBOOL           dodry       = TRUE;   // do dry dynamics?
          GBOOL           docoriolis  = FALSE;  // use Coriolis force?
          GBOOL           dograv      = TRUE;   // use gravitational force?
          GBOOL           dofallout   = FALSE;  // allow precip fallout?
          GBOOL           bconserved  = FALSE;  // use conserved form?
          GBOOL           bforced     = FALSE;  // use forcing?
          GBOOL           usemomden   = TRUE;   // use momentum density form?
          GBOOL           usebase     = TRUE;   // use hydrostatic base state?
          GBOOL           variabledt  = FALSE;  // use variable timestep?
          GBOOL           bvarvterm   = FALSE;  // time-variable term vel?
          GINT            nstate      = GDIM+2; // no. vars in state vec
          GINT            nsolve      = GDIM+2; // no. vars to solve for
          GINT            nlsector    = 0;      // no. vars in liq-sector
          GINT            nisector    = 0;      // no. vars in ice-sector
          GINT            nbase       = 2;      // no. vars in base state
          GINT            ntmp        = 8;
          GINT            itorder     = 2;
          GINT            inorder     = 2;
          GStepperType    isteptype   = GSTEPPER_EXRK;
          GFTYPE          courant     = 0.5;    // Courant factor
          GFTYPE          nu          = 0.0;    // viscosity constant
          GTVector<GINT>  iforced;              // state comps to foce
          GTVector<Ftype> omega;                // rotation rate vector
          GString         ssteptype;            // stepping method
        };

        GMConv() = delete; 
        GMConv(Grid &grid, GMConv<TypePack>::Traits &traits, State &tmp);
       ~GMConv();
        GMConv(const GMConv &bu) = default;
        GMConv &operator=(const GMConv &bu) = default;

        StateComp           &get_nu() { return nu_; };                       // Set nu/viscosity

        void                set_steptop_callback(
                            std::function<void(const Time &t, State &u, 
                                               const Time &dt)> callback) 
                             { steptop_callback_ = callback; bsteptop_ = TRUE;}
                                            

protected:
        void                step_impl(const Time &t, State &uin, State &uf, State &ub, 
                                      const Time &dt);                    // Take a step
        void                step_impl(const Time &t, const State &uin, State &uf, State &ub,
                                      const Time &dt, State &uout);       // Take a step
        GBOOL               has_dt_impl() const {return traits_.variabledt;}  
                                                                          // Has dynamic dt?
        void                dt_impl(const Time &t, State &u, Time &dt);   // Get dt
        void                apply_bc_impl(const Time &t, State &u, 
                                          State &ub);                     // Apply bdy conditions
private:

        void                init();                                       // initialize 
        GINT                req_tmp_size();                               // required tmp size
        void                dudt_impl  (const Time &t, const State &u, const State &uf, const State &ub,
                                        const Time &dt, Derivative &dudt);
        void                step_exrk  (const Time &t, State &uin, State &uf, State &ub,
                                        const Time &dt, State &uout);
        void                step_multistep(const Time &t, State &uin, State &uf, State &ub,
                                           const Time &dt);
        void                cycle_keep   (const State &u);
inline  void                compute_cv   (const State &u, State &utmp, StateComp &cv);
inline  void                compute_qd   (const State &u, State &utmp, StateComp &qd);
inline  void                compute_temp (const State &u, State &utmp, StateComp &t );
inline  void                compute_p    (const State &u, State &utmp, StateComp &p );
inline  void                compute_ptemp(const State &u, State &utmp, StateComp &temp, StateComp &p );
inline  void                compute_falloutsrc
                                         (StateComp &g, State &qi, State &v, GINT jexcl, State &utmp, StateComp &r );
inline  void                compute_div  (StateComp &q, State &v, State &utmp, StateComp &div );
inline  void                compute_v    (const State &u, State &utmp, State &v);
inline  void                compute_vterm(StateComp &tv, State &W);
inline  void                assign_helpers(const State &u, const State &uf);
inline  void                compute_pe   (StateComp &rhoT, State &qi, State &tvi, State &utmp, StateComp      &r);
 

        GBOOL               bforced_;       // use forcing vectors
        GBOOL               bsteptop_;      // is there a top-of-step callback?
        GBOOL               bvterm_;        // teminal vel. computed?
        GINT                nevolve_;       // num StateComp's evolved
        GINT                nhydro_;        // num hydrometeors
        GINT                nmoist_;        // number of moist components
        GStepperType        isteptype_;     // stepper type
        GTVector<GFTYPE>    tcoeffs_;       // coeffs for time deriv
        GTVector<GFTYPE>    acoeffs_;       // coeffs for NL adv term
        GTVector<GFTYPE>    dthist_;        // coeffs for NL adv term
        State               uevolve_;       // helper array to specify evolved sstate components
        State               utmp_;
        State               uold_;          // helper arrays set from utmp
        State               urhstmp_;       // helper arrays set from utmp
        State               uoptmp_;        // helper arrays set from utmp
        State               urktmp_;        // helper arrays set from utmp
        State               fk_;            // kinetic forcing vector
        State               qi_;            // full mass fraction vector
        State               qice_;          // ice mass fraction vector
        State               qliq_;          // liquid mass fraction vector
        State               tvi_;           // term vel. vector for all qi
        State               tvice_;         // term vel. vector for all qice
        State               tvliq_;         // term vel. vector for all qliq
        State               fv_;            // state velocity forcing components
        State               s_;             // state momentum components
        State               v_;             // state velocity components
        State               W_;             // terminal velocity components
        GTVector<State>     ukeep_;         // state at prev. time levels
        GTVector<GString>
                            valid_types_;   // valid stepping methods supported
        GTVector<GFTYPE>    nu_   ;         // dissipoation
        GTVector<GFTYPE>    dxmin_ ;        // element face mins
        GTVector<GFTYPE>    maxbyelem_ ;    // element-based maxima for dt
        GGrid              *grid_;          // GGrid object
        GExRKStepper<GFTYPE>
                           *gexrk_;         // ExRK stepper, if needed
        GMass              *gmass_;         // mass op
        GMass              *gimass_;        // inverse mass op
        GAdvect            *gadvect_;       // advection op
        GHelmholtz         *ghelm_;         // Helmholz and Laplacian op
        GpdV<TypePack>     *gpdv_;          // pdV op
//      GFlux              *gflux_;         // flux op
        GC_COMM             comm_;          // communicator
        GGFX<GFTYPE>       *ggfx_;          // gather-scatter operator
        GMConv<TypePack>::Traits 
                            traits_;        // traits structure

        std::function<void(const Time &t, State &u, const Time &dt)>
                           steptop_callback_;

        // Starting indices for each sector:
        GINT               MOMENTUM;
        GINT               ENERGY;
        GINT               DENSITY;
        GINT               VAPOR;
        GINT               LIQMASS;
        GINT               ICEMASS;
        GINT               PRESCRIBED;
        GINT               BASESTATE;
        GINT               LIQTERMV;
        GINT               ICETERMV;


};

#include "gmconv.ipp"

#endif
