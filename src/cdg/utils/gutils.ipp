//==================================================================================
// Module       : gutils.ipp
// Date         : 1/31/19 (DLR)
// Description  : GeoFLOW utilities namespace
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

namespace geoflow
{


//**********************************************************************************
//**********************************************************************************
// METHOD : smooth
// DESC   :
//          
// DESC   : Computes a weighted average.
//              Calculates:
//                u <- DSS(M_L u) / DSS(M_L),
//          where u is the field expressed in local format;
//          M_L is the local mass matrix (unassembled);
//          DSS is the application of Q Q^T, or the direct stiffness operator.
// ARGS   : 
//          grid : GGrid object
//          tmp  : tmp space 
//          op   : GGFX_OP 
//          u    : Locally expressed field to smooth
// RETURNS: none.
//************************************************************************************
template<typename T>
void smooth(GGrid &grid, GGFX_OP op, GTVector<T> &tmp, GTVector<T> &u)
{
  GGFX *ggfx = &grid->get_ggfx();
  tmp = u;
 
  u.pointProd(grid_->massop()->data());
  tmp = grid_->imassop()->data();
  ggfx_->doOp(tmp, op);  // DSS(mass_local)

  u.pointProd(tmp);      // M_assembled u M_L

} // end, smooth method


} // end, namespace

