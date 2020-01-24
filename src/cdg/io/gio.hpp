//==================================================================================
// Module       : gio.hpp
// Date         : 1/20/20 (DLR)
// Description  : GIO object encapsulating methods for POSIX and collective IO
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : IOBase.
//==================================================================================
#if !defined(_GIO_HPP)
#define _GIO_HPP

#include "gtvector.hpp"
#include "ggrid.hpp"
#include "pdeint/io_base.hpp"
#include "tbox/property_tree.hpp"
#include "tbox/mpixx.hpp"

using namespace geoflow::pdeint;
using namespace std;


template<typename IOType>
class GIO : public IOBase<IOType>
{

public:
        enum GIOType       {GIO_POSIX=0, GIO_COLL}; // POSIX or collective
        using Types       = IOType;
        using IOBaseType  = IOBase<IOType>;
        using IOBasePtr   = std::shared_ptr<IOBaseType>;
        using State       = typename Types::State;
        using Grid        = typename Types::Grid;
        using Value       = typename Types::Value;
        using Time        = typename Types::Time;
        using Size        = typename Types::Size;
        using StateInfo   = typename Types::StateInfo; // May contain time, time index, IO dir

        struct Traits {
          GIOType     io_type = GIO_COLL; // default to collective IO
          GINT        ivers   = 0;        // IO version tag
          GBOOL       multivar= false;    // multiple vars in file (only of COLL types)?
          GINT        wtime   = 6;        // time-field width
          GINT        wtask   = 5;        // task-field width (only for POSIX types)
          GINT        wfile   = 2048;     // file name max
          GINT        dim     = GDIM;     // problem dimension
        }; 

        struct GIOStateTraits {
          GINT        gtype   = 0;        // check src/cdg/include/gtypes.h
          GINT        sttype  = 1;        // state type index (grid=0 or state=1)
          GSIZET      index   = 0;        // time index
          GSIZET      nelems  = 0;        // num elems
          GSIZET      cycle   = 0;        // continuous time cycle
          GFTYPE      time    = 0.0;      // state time
          std::vector<GINT> 
                      indirect;           // which elements of state to do IO on
          std::vector<GString> 
                      svars;              // names of state members
          GTMatrix<GINT>
                      porder;             // if ivers=0, is 1 X GDIM; else nelems X GDIM;
          GString     idir    = ".";      // input directory
          GString     odir    = ".";      // output directory
        }; 

        static_assert(std::is_same<Value,GFTYPE>::value,
               "Value is of incorrect type");
        static_assert(std::is_same<State,GTVector<GTVector<GFTYPE>*>>::value,
               "State is of incorrect type");
        static_assert(std::is_same<StateInfo,GStateIOTraits>::value,
               "StateInfo is of incorrect type");

                           GIO() = delete;
                           GIO(Grid &grid, Traits &traits, GC_COMM comm);
                          ~GIO() = default;
                           GIO(const GIO &a) = default;
                           GIO &operator=(const GIO &bu) = default;

        void               write_state_impl(std::string filename, StateInfo &info, State &u);
        void               read_state_impl(std::string filename, StateInfo &info, State &u);

private:
// Private methods:
        void               init(const Time &t, const State &u);
        void               write_state_posix(StateInfo &info, const State  &u);
        void               read_state_posix (StateInfo &info,       State  &u);
        void               write_grid_posix (StateInfo &info);
        GSIZET             write_posix(GString filename, StateInfo &info, const GTVector<Value> &u);
        GSIZET             read_posix (GString filename, StateInfo &info,       GTVector<Value> &u);
        GSIZET             read_header(GString filename, StateInfo &info, Traits &traits);
        void               resize(GINT n);


// Private data:
        GBOOL              bInit_;      // is initialized?
        GINT               nfname_;
        GC_COMM            comm_;
        GString            default_state_name_pref_;
        GString            default_grid_name_pref_;
        std::stringstream  svarname_;
        std::stringstream  sgridname_;
        GString            fname_;
        char              *cfname_;
        Traits            *traits_;
        Grid              *grid_;
        GTVector<GINT>     state_index_;// list of state indices to print
        GTVector<GString>  state_names_;// list of names of state files
        GTVector<GString>  grid_names_ ;// list of names of grid comp files

};

#include "gio.ipp"

#endif

