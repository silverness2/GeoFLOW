//==================================================================================
// Module       : ggfx.ipp
// Date         : 5/9/18 (DLR)
// Description  : Template functions implementations for GGFX
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

#include <math.h>
#include <type_traits>
#include <assert.h>
#include <limits>
#include "ggfx.hpp"
#include "gcomm.hpp"


//************************************************************************************
//************************************************************************************
// METHOD : DoOp (1)
// DESC   : Actually carries out the operation, op, on the elements
//          of field, u, using the intialization from the call to Init.
// ARGS   : u    : input field whose elements are to be combined. Must
//                 be same size as glob_index in Init call.
//          nu   : size of buffer, u
//          op   : operation to be performed in the combination
// RETURNS: TRUE on success; else FALSE, if invalid operation requested, or if 
//          there is an error in the gather/scatter combination.
//************************************************************************************
template<typename T> 
GBOOL GGFX::DoOp(T *&u, GSIZET nu,  GGFX_OP op)
{
  assert( std::is_arithmetic<T>::value && "Illegal template type");

  GINT      irank;
  GINT      i, j;
  GBOOL     bret=TRUE;
  GString   serr = "GGFX::DoOp(1): ";

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_, serr << "Doing first LocalGS...");
#endif
  
  // For each global index row in iOpL2LIndices, gather and
  // scatter data locally: get one operation for each of the
  // shared nodes; 
  bret = LocalGS(u, nu, iOpL2LIndices_, nOpL2LIndices_, op);
  if ( !bret ) {
    std::cout << serr << "LocalGS (1) failed " << std::endl;
    exit(1);
  }
  
#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "First LocalGS done.");
#endif

#if defined(GGFX_TRACE_OUTPUT)
  GTVector<T> utmp(u,nu);
  GPP(comm_,serr << "After LocalGS: u=" << utmp);
  GPP(comm_,serr << "iL2RIndices=" << iOpL2RIndices_);
#endif

  // If no other tasks, return:
  if ( nprocs_ == 1 ) bret = TRUE;

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "Doing Data_Exchange...");
#endif

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "Doing Data_Exchange");
#endif
  // Perform a exchange of field data:
  bret = Data_Exchange(u, nu);
  if ( !bret ) {
    std::cout << serr << "Data_Exchange failed " << std::endl;
    exit(1);
  }
#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "Data_Exchange done.");
#endif

  bret = LocalGS(u, nu, iOpR2LIndices_, nOpR2LIndices_, op, recvBuff_.data().data());
  if ( !bret ) {
    std::cout << serr << "LocalGS (2) failed " << std::endl;
    exit(1);
  }

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "done.");
#endif
  return bret;

} // end of method DoOp(1)


//************************************************************************************
//************************************************************************************
// METHOD : DoOp (1)
// DESC   : Actually carries out the operation, op, on the elements
//          of field, u, using the intialization from the call to Init.
// ARGS   : u    : input field whose elements are to be combined. Must
//                 be same size as glob_index in Init call.
//          op   : operation to be performed in the combination
// RETURNS: TRUE on success; else FALSE, if invalid operation requested, or if 
//          there is an error in the gather/scatter combination.
//************************************************************************************
template<typename T>  
GBOOL GGFX::DoOp(GTVector<T> &u, GGFX_OP op) 
{
   GSIZET   nu = u.size();
   T       *uu = u.data();
   return DoOp(uu, nu,  op);

} // end, method DoOp (1)


//************************************************************************************
//************************************************************************************
// METHOD : DoOp (2)
// DESC   : Actually carries out the operation, op, on the elements
//          of field, u, using the intialization from the call to Init.
// ARGS   : u    : input field whose elements are to be combined. Can be 
//                 _any_ size, but the number of elems use din iind array 
//                 _must_ be the same size as the number in glob_index in Init call. 
//          nu   : size of buffer, u
//          iind : indirection array into u. Must contain indices that
//                 corresp to those in glob_index in Init call.
//          nind : number of indirection indices in iind; must be the same
//                 as the number of glob_index in Init call.
//          op   : operation to be performed in the combination
// RETURNS: TRUE on success; else FALSE, if invalid operation requested, or if 
//          there is an error in the gather/scatter combination.
//************************************************************************************
template<typename T> GBOOL GGFX::DoOp(T *&u, GSIZET nu,  GSIZET *iind, GSIZET nind, GGFX_OP op)
{
  assert( std::is_arithmetic<T>::value && "Illegal template type");

  GINT      irank;
  GINT      i, j;
  GBOOL     bret=TRUE;
  GString   serr = "GGFX::DoOp(2): ";

#if defined(GGFX_TRACE_OUTPUT)
  std::cout << serr << "Doing first LocalGS..." << std::endl;
#endif

  // For each global index row in iOpL2LIndices, gather and
  // scatter data locally: get one operation for each of the
  // shared nodes; 
  bret = LocalGS(u, nu, iind, nind, iOpL2LIndices_, nOpL2LIndices_, op);
  if ( !bret ) {
    std::cout << serr << "LocalGS (1) failed " << std::endl;
    exit(1);
  }
  
#if defined(GGFX_TRACE_OUTPUT)
  std::cout << serr << "First LocalGS done ." << std::endl;
#endif

  // If no other tasks, return:
  if ( nprocs_ == 1 ) bret = TRUE;


#if defined(GGFX_TRACE_OUTPUT)
  std::cout << serr << "Doing Data_Exchange..." << std::endl;
#endif

  // Perform a exchange of field data:
  bret = Data_Exchange(u, nu, iind, nind);
  if ( !bret ) {
    std::cout << serr << "Data_Exchange failed " << std::endl;
    exit(1);
  }
#if defined(GGFX_TRACE_OUTPUT)
  std::cout << serr << "Data_Exchange done ." << std::endl;
#endif


#if defined(GGFX_TRACE_OUTPUT)
  std::cout << serr << "Doing second LocalGS..." << std::endl;
#endif
  bret = LocalGS(u, nu, iind, nind, iOpL2RIndices_, nOpL2RIndices_, op, recvBuff_.data().data());
  if ( !bret ) {
    std::cout << serr << "LocalGS (2) failed " << std::endl;
    exit(1);
  }
#if defined(GGFX_TRACE_OUTPUT)
  std::cout << serr << "Second LocalGS done." << std::endl;
#endif

#if defined(GGFX_TRACE_OUTPUT)
  std::cout << serr << "Doing final assignment..." << std::endl;
#endif

#if defined(GGFX_TRACE_OUTPUT)
  std::cout << serr << "done." << std::endl;
#endif
  return bret;

} // end of method DoOp(2)


//************************************************************************************
//************************************************************************************
// METHOD : DoOp (2)
// DESC   : Actually carries out the operation, op, on the elements
//          of field, u, using the intialization from the call to Init.
// ARGS   : u    : input field whose elements are to be combined. Can be 
//                 _any_ size, but the number of elems use din iind array 
//                 _must_ be the same size as the number in glob_index in Init call. 
//          iind : indirection array into u. Must contain indices that
//                 corresp to those in glob_index in Init call.
//          op   : operation to be performed in the combination
// RETURNS: TRUE on success; else FALSE, if invalid operation requested, or if 
//          there is an error in the gather/scatter combination.
//************************************************************************************
template<typename T>  
GBOOL GGFX::DoOp (GTVector<T> &u,  GSZBuffer &ind, GGFX_OP op)
{
   GSIZET   nind = ind.size();
   GSIZET   iind = ind.data();
   GSIZET   nu = u.size();
   T       *uu = u.data();
   return DoOp(uu, nu, iind, nind, op);

} // end, method DoOp (2)


//************************************************************************************
//*************************************************************************************
// METHOD : LocalGS (1)
// DESC   : Performs GGFX_OP (op) on the array of field values, u
//          and scatters them to the same dof
// ARGS   : u      : input field whose elements are to be combined. May
//                   be larger than number of elements in glob_index array
//                   used in call to Init, as long as indirection indices,
//                   iind, are being used
//          nu     : number of elements in u
//          ilocal : Matrix of indirection indices into qu array. Each column represents
//                   a shared node, and each row represents the other local indices 
//                   that point to the same shared node, so that the values at the same
//                   global index can be summed, multiplied, max'd or min'd.
//          nlocal : for each shared node (column in ilocal), local indices to it
//          op     : operation to be performed in the combination
//          qop    : operand; used if non_NULL; must have same no elements as columns in
//                   ilocal.
// RETURNS: TRUE on success; else FALSE
//************************************************************************************
template<typename T> 
GBOOL GGFX::LocalGS(T *&qu, GSIZET nu, GSZMatrix &ilocal, GSZBuffer &nlocal, GGFX_OP op, GDOUBLE *qop)
{

  assert( std::is_arithmetic<T>::value && "Illegal template type");

  GString serr = "GGFX::LocalGS (1): ";
  T       res;
  GLLONG  i, j;


  // Perform GGFX_OP on the nodes shared by this proc:
  switch(op) {
    case GGFX_OP_SUM:
      if ( qop == NULLPTR ) {
        for ( j=0; j<ilocal.size(2); j++ ) {
          res = 0.0;
          for ( i=0; i<nlocal[j]; i++ ) { // do gather
             res += qu[ilocal(i,j)];
          }
          for ( i=0; i<nlocal[j]; i++ ) { // do scatter
             qu[ilocal(i,j)] = res;
          }
        }
      }
      else {
        for ( j=0; j<ilocal.size(2); j++ ) {
          for ( i=0; i<nlocal[j]; i++ ) { // do scatter
             qu[ilocal(i,j)] += qop[j];
          }
        }
      }
      break;
    case GGFX_OP_PROD:
      if ( qop == NULLPTR ) {
        for ( j=0; j<ilocal.size(2); j++ ) {
          res = 1.0;
          for ( i=0; i<nlocal[j]; i++ ) {
             res *= qu[ilocal(i,j)];
          }
          for ( i=0; i<nlocal[j]; i++ ) {
             qu[ilocal(i,j)] = res;
          }
        }
      }
      else {
        for ( j=0; j<ilocal.size(2); j++ ) {
          for ( i=0; i<nlocal[j]; i++ ) {
             qu[ilocal(i,j)] *= qop[j];
          }
        }
      }
      break;
    case GGFX_OP_MAX:
      if ( qop == NULLPTR ) {
        for ( j=0; j<ilocal.size(2); j++ ) {
          res = std::numeric_limits<T>::min();
          for ( i=0; i<nlocal[j]; i++ ) {
             res = MAX(res,qu[ilocal(i,j)]);
          }
          for ( i=0; i<nlocal[j]; i++ ) {
             qu[ilocal(i,j)] = res;
          }
        }
      }
      else {
        for ( j=0; j<ilocal.size(2); j++ ) {
          for ( i=0; i<nlocal[j]; i++ ) {
             qu[ilocal(i,j)] = MIN(qu[ilocal(i,j)],qop[j]);
          }
        }
      }
      break;
    case GGFX_OP_MIN:
      if ( qop == NULLPTR ) {
        for ( j=0; j<ilocal.size(2); j++ ) {
          res = std::numeric_limits<T>::max();
          for ( i=0; i<nlocal[j]; i++ ) {
             res = MIN(res,qu[ilocal(i,j)]);
          }
          for ( i=0; i<nlocal[j]; i++ ) {
             qu[ilocal(i,j)] = res;
          }
        }
      }
      else {
        for ( j=0; j<ilocal.size(2); j++ ) {
          res = std::numeric_limits<T>::max();
          for ( i=0; i<nlocal[j]; i++ ) {
             qu[ilocal(i,j)] = MIN(qu[ilocal(i,j)],qop[j]);
          }
        }
      }
      break;
    default:
      return FALSE;
  }

  return TRUE;

} // end of method LocalGS (1)


//************************************************************************************
//*************************************************************************************
// METHOD : LocalGS (2)
// DESC   : Performs GGFX_OP (op) on the array of field values, u,
//          and scatters them to the same dof
// ARGS   : u      : input field whose elements are to be combined. May
//                   be larger than number of elements in glob_index array
//                   used in call to Init, as long as indirection indices,
//                   iind, are being used
//          nu     : number of elements in u
//          iind   : indirection array into u. Must contain indices that
//                   corresp to those in glob_index in Init call.
//          nind   : number of indirection indices in iind; must be the same
//                   as the number of glob_index in Init call.
//          ilocal : Matrix of indirection indices into qu array. Each row represents
//                   a shared node, and each column represents the other local indices 
//                   that point to the same shared node, so that the values at the same
//                   global index can be summed, multiplied, max'd or min'd.
//          nlocal : for each shared node (row in ilocal), local indices to it
//          op     : operation to be performed in the combination
//          qop    : operand; used if non_NULL; must have same no elements as columns in
//                   ilocal.
// RETURNS: TRUE on success; else FALSE
//************************************************************************************
template<typename T> 
GBOOL GGFX::LocalGS(T *&qu, GSIZET nu, GSIZET *&iind, GSIZET nind, 
                    GSZMatrix &ilocal, GSZBuffer &nlocal,  GGFX_OP op, GDOUBLE *qop)
{

  assert( std::is_arithmetic<T>::value && "Illegal template type");

  GString serr = "GGFX::LocalGS (2): ";
  T       res;
  GLLONG  i, j;


  // Perform GGFX_OP on the nodes shared by this proc:
  switch(op) {
    case GGFX_OP_SUM:
      if ( qop == NULLPTR ) {
        for ( j=0; j<ilocal.size(2); j++ ) {
          res = 0.0;
          for ( i=0; i<nlocal[j] && ilocal(i,j)<nind; i++ ) {
             res += qu[iind[ilocal(i,j)]];
          }
          for ( i=0; i<nlocal[j]; i++ ) {
             qu[iind[ilocal(i,j)]] = res;
          }
        }
      }
      else {
        for ( j=0; j<ilocal.size(2); j++ ) {
          for ( i=0; i<nlocal[j]; i++ ) {
             qu[iind[ilocal(i,j)]] += qop[j];
          }
        }
      }
      break;
    case GGFX_OP_PROD:
      if ( qop == NULLPTR ) {
        for ( j=0; j<ilocal.size(2); j++ ) {
          res = 0.0;
          for ( i=0; i<nlocal[j] && ilocal(i,j)<nind; i++ ) {
             res *= qu[iind[ilocal(i,j)]];
          }
          for ( i=0; i<nlocal[j]; i++ ) {
             qu[iind[ilocal(i,j)]] = res;
          }
        }
      }
      else {
        for ( j=0; j<ilocal.size(2); j++ ) {
          for ( i=0; i<nlocal[j]; i++ ) {
             qu[iind[ilocal(i,j)]] *= qop[j];
          }
        }
      }
      break;
    case GGFX_OP_MAX:
      if ( qop == NULLPTR ) {
        for ( j=0; j<ilocal.size(2); j++ ) {
          res = std::numeric_limits<T>::min();
          for ( i=0; i<nlocal[j] && ilocal(i,j)<nind; i++ ) {
             res = MAX(res,qu[iind[ilocal(i,j)]]);
          }
          for ( i=0; i<nlocal[j]; i++ ) {
             qu[iind[ilocal(i,j)]] = res;
          }
        }
      }
      else {
        for ( j=0; j<ilocal.size(2); j++ ) {
          for ( i=0; i<nlocal[j]; i++ ) {
             qu[iind[ilocal(i,j)]] = MAX(qu[iind[ilocal(i,j)]],qop[j]);
          }
        }
      }
      break;
    case GGFX_OP_MIN:
      if ( qop == NULLPTR ) {
        for ( j=0; j<ilocal.size(2); j++ ) {
          res = std::numeric_limits<T>::max();
          for ( i=0; i<nlocal[j]&& ilocal(i,j)<nind; i++ ) {
             res = MIN(res,qu[iind[ilocal(i,j)]]);
          }
          for ( i=0; i<nlocal[j]; i++ ) {
             qu[iind[ilocal(i,j)]] = res;
          }
        }
      }
      else {
        for ( j=0; j<ilocal.size(2); j++ ) {
          for ( i=0; i<nlocal[j]; i++ ) {
             qu[iind[ilocal(i,j)]] = MAX(qu[iind[ilocal(i,j)]],qop[j]);
          }
        }
      }
      break;
    default:
      return FALSE;
  }

  return TRUE;

} // end of method LocalGS (2)


//************************************************************************************
//************************************************************************************
// METHOD : Data_Exchange (1)
// DESC   : Perform data exchange between procs by using asynchronous
//              recvs...
// ARGS   :
//          u    : field data 
//          nu   : field data received. Allocated in Init.
// RETURNS: TRUE on success; else FALSE
//************************************************************************************
template<typename T> 
GBOOL GGFX::Data_Exchange(T *&u, GSIZET nu)
{
  assert( std::is_arithmetic<T>::value && "Illegal template type");

  GString       serr = "GGFX::Data_Exchange (1): ";
  GCommDatatype dtype=T2GCDatatype<GDOUBLE>();
  GINT          i, j;
  GBOOL         bret=TRUE;

  sendBuff_.set(-1);
  recvBuff_.set(-1);

  // Fill send buffer with data:
  for ( j=0; j<iOpL2RTasks_.size() && bret; j++ ) {
 // if ( iOpL2RTasks_[j] == rank_ ) continue;
    for ( i=0; i<nOpL2RIndices_[j]; i++ ) {
      sendBuff_(i,j) = u[iOpL2RIndices_(i,j)];
    }
  }

#if defined(GGFX_TRACE_OUTPUT)
  GIBuffer rtasks(iOpL2RTasks_);
  GIBuffer isend(iOpL2RTasks_.size());
  GPP(comm_,serr << "iOpL2RTasks=" << iOpL2RTasks_);
  GPP(comm_,serr << "sendBuff=" << sendBuff_);
#endif

  bret = GComm::ASendRecv(recvBuff_.data().data(),iOpL2RTasks_.size(),NULLPTR,recvBuff_.size(1),dtype,iOpL2RTasks_.data(), TRUE, 
                          sendBuff_.data().data(),iOpL2RTasks_.size(),NULLPTR,sendBuff_.size(1),dtype,iOpL2RTasks_.data(), comm_);

#if defined(GGFX_TRACE_OUTPUT)
  GPP(comm_,serr << "(2)sendBuff=" << sendBuff_);
  GPP(comm_,serr << "(2)recvBuff=" << recvBuff_);
#endif

  return bret;
} // end of Data_Exchange (1)


//************************************************************************************
//************************************************************************************
// METHOD : Data_Exchange (2)
// DESC   : Perform data exchange between procs by using asynchronous
//              recvs...
// ARGS   :
//          u    : field data 
//          nu   : field data received. Allocated in Init.
//          iind : indirection array into u. Must contain indices that
//                 corresp to those in glob_index in Init call.
//          nind : number of indirection indices in iind; must be the same
//                 as the number of glob_index in Init call.
// RETURNS: TRUE on success; else FALSE
//************************************************************************************
template<typename T> 
GBOOL GGFX::Data_Exchange(T *&u, GSIZET nu, GSIZET *&iind, GSIZET nind)
{
  assert( std::is_arithmetic<T>::value && "Illegal template type");

  GString       serr = "GGFX::Data_Exchange (2): ";
  GCommDatatype dtype=T2GCDatatype<GDOUBLE>();
  GINT           i, j;
  GBOOL          bret=TRUE;


  // Fill send buffer with data:
  for ( j=0; j<iOpL2RTasks_.size() && bret; j++ ) {
//  if ( iOpL2RTasks_[j] == rank_ ) continue;
    for ( i=0; j<nOpL2RIndices_[j] && iOpL2RIndices_(i,j)<nind; i++ ) {
      sendBuff_(i,j) = u[iind[iOpL2RIndices_(i,j)]];
    }
  }

recvBuff_.set(-1);

  bret = GComm::ASendRecv(recvBuff_.data().data(),(GINT)recvBuff_.size(2),NULLPTR,(GINT)recvBuff_.size(1),dtype,iOpL2RTasks_.data(), TRUE, 
                          sendBuff_.data().data(),(GINT)sendBuff_.size(2),NULLPTR,(GINT)sendBuff_.size(1),dtype,iOpL2RTasks_.data(), comm_);


  return bret;
} // end of Data_Exchange (2)

