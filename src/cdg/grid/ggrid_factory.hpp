//==================================================================================
// Module       : ggrid_factory
// Date         : 2/1/19 (DLR)
// Description  : GeoFLOW grid factory object. 
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================
#if !defined(_GGRID_FACTORY_HPP)
#define _GGRID_FACTORY_HPP 

#include "tbox/property_tree.hpp"
#include "gcomm.hpp"
#include "gtvector.hpp"
#include "ggrid.hpp"
#include "pdeint/io_base.hpp"
#include "pdeint/observer_base.hpp"
#include "ggrid_icos.hpp"
#include "ggrid_box.hpp"


typedef GFTYPE                      Time;
typedef GTVector<GTVector<GFTYPE>*> State;

template<typename TypePack>
class GGridFactory
{
  public:

        using StateInfo   = typename TypePack::StateInfo;
        using IOBaseType  = IOBase<TypePack>;
        using IOBasePtr   = std::shared_ptr<IOBaseType>;
        using Grid        = typename TypePack::Grid;

	static GGrid *build(const geoflow::tbox::PropertyTree& ptree, GTVector<GNBasis<GCTYPE,GFTYPE>*> gbasis, IOBasePtr pIO, typename ObserverBase<TypePack>::Traits &obstraits, GC_COMM &comm);


  private:
        static void   read_grid(const geoflow::tbox::PropertyTree& ptree, GTMatrix<GINT> &p, GTVector<GTVector<GFTYPE>> &xnodes, IOBasePtr pIO, typename ObserverBase<TypePack>::Traits &obstraits, GC_COMM &comm);


}; // class GGridFactory

#include "ggrid_factory.ipp"

#endif 
