//==================================================================================
// Module       : gio.hpp
// Date         : 1/20/20 (DLR)
// Description  : GIO object encapsulating methods for POSIX and collective IO
// Copyright    : Copyright 2020. Colorado State University. All rights reserved.
// Derived From : IOBase.
//==================================================================================

//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method (1)
// DESC   : Instantiate with Traits
// ARGS   : grid  : grid object
//          traits: Traits sturcture
//          comm  : communicator
//**********************************************************************************
template<typename IOType>
GIO<IOType>::GIO(Grid &grid,  Traits &traits, GC_COMM comm):
IOBase<IOType>(grid, traits),
nstates_                       (0),
comm_                       (comm),
cfname_                  (NULLPTR),
nfname_                        (0)
{ 
#if !defined(_G_USE_MPI)
  assert(this->traits_->io_type == IOBase<IOType>::GIO_COLL && "Collective IO only allowed if MPI is used");
#endif
  GSIZET           ndof  =this->grid_->ndof();

  extent_.resize(GComm::WorldSize(comm_));

  // Get extents of a single state component on each task:
  GCom::Allgather(&ndof, 1, T2GCDatatype<GSIZET>(), extent_.data(), 1, T2GCDatatype<GSIZET>(), comm_);

} // end of constructor (1) method


//**********************************************************************************
//**********************************************************************************
// METHOD : Destructor method 
// DESC   :
// ARGS   : none.
//**********************************************************************************
template<typename IOType>
GIO<IOType>::~GIO()
{ 
#if defined(_G_USE_MPI)
  MPI_Type_free(&mpi_state_type_);
#endif
} // end of destructor method


//**********************************************************************************
//**********************************************************************************
// METHOD     : init
// DESCRIPTION: Initialize object
// ARGUMENTS  : t  : state time
//              u  : state variable
// RETURNS    : none.
//**********************************************************************************
template<typename IOType>
void GIO<IOType>::init()
{

  if ( this->traits_.io_type != IOBase<IOType>::GIO_COLL ) {
    return; // nothing to do
  }

  GSIZET           ndof  =this->grid_->ndof();
  GSIZET           nelems=this->grid_->nelems();

  GINT myrank = GComm::WorldRank(comm_);
  if ( myrank == 0 ) {
    state_disp_ = 0; // in bytes
    state_extent_ = extent[myrank]*sizeof(GFTYPE); // in bytes
  }
  else {
    state_disp_   = extent.sum(0,myrank)*sizeof(GFTYPE); // in bytes
    state_extent_ = extent[myrank]*sizeof(GFTYPE); // in bytes
  }
  nbgdof_ = extent.sum()*sizeof(GFTYPE); // no. bytes of state

#if defined(_G_USE_MPI)
  MPI::Aint lowerbnd=0;
  MPI_Type_create_resized(T2GCDatatype<GFTYPE>(), lowerbnd, state_extent_, &mpi_state_type_);
  MPI_Type_commit(&mpi_state_type_);

#endif
  bInit_ = TRUE;


} // end of method init



//**********************************************************************************
//**********************************************************************************
// METHOD : write_state_impl
// DESC   : Do GIO POSIX or collective output of state. 
//          If there aren't enough svar specified, or if filename isn't specified,
//          when needed, then it's an error.
//          Note: if info.sttype > 0, then we assume we are printing a 
//          grid, and the filename is created without the time index tag.
// ARGS   : filename: used if info.multivar > 0 to specify file name for
//                    all state variables. This works only if we GIOType is GIO_COLL.
//                    If info.multivar ==0, individual filename refixes are provided 
//                    in info.svar
//          info    : StateInfo structure
//          u       : state
// RETURNS: none
//**********************************************************************************
template<typename IOType>
void GIO<IOType>::write_state_impl(std::string filename, StateInfo &info, const State &u)
{
    GString        serr = "write_state_impl: ";
    GINT           myrank = GComm::WorldRank(comm_);
    GSIZET         nb, nc, nd;
    GTVector<GTVector<GFTYPE>>
                  *xnodes = &(this->grid_->xNodes());
    State         ostate(1);
    GElemList     *elems  = &(this->grid_->elems());
    std::stringstream 
                   cformat, scformat, spformat;

    if ( !bInit_ ) init();

    // Required number of coord vectors:
    nc = this->grid_->gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;
    if ( info.sttype !=0 ) { // is a grid file 
      assert(nc ==  xnodes->size()
          && "Incompatible grid or coord  dimensions");
    }
    
    this->traits_.dim  = GDIM;

    if ( this->traits_.io_type == IOBase<IOType>::GIO_COLL ) {
      info.nelems = this->grid_->ngelems(); // total no. elems among all tasks
    }
    else {
      info.nelems = this->grid_->nelems(); // local nelems for this task
    }
    info.gtype    = this->grid_->gtype();


    // Build format strings:
    resize(this->traits_.wfile);
    scformat.str(""); spformat.clear();
    spformat.str(""); spformat.clear();
    if ( info.sttype ==0 ) { // is a physical state
      spformat   << "%s/%s.%0" << this->traits_.wtime << "d.%0" << this->traits_.wtask << "d.out";
      scformat   << "%s/%s.%0" << this->traits_.wtask << "d.out";
    }
    else {                   // is a grid file
      spformat   << "%s/%s.%0" << this->traits_.wtask << "d.out";
      scformat   << "%s/%s.%0" << this->traits_.wtask << "d.out";
    }

    // Set porder vector depending on version:
    info.porder.resize(this->traits_.ivers == 0 ? 1 : info.nelems, GDIM);
    for ( auto i=0; i<info.porder.size(1); i++ )  { // for each element
      for ( auto j=0; j<info.porder.size(2); j++ ) info.porder(i,j) = (*elems)[i]->order(j);
    }

    if ( !this->traits_.multivar ) { // one state comp per file
      assert(info.svars.size() >= u.size());
      // Cycle over all fields, and write:
      for ( auto j=0; j<u.size(); j++ ) {
        assert(u[j]->size() > 0 && "Invalid state component");
        svarname_.str(""); svarname_.clear();
        assert(info.svars[j].length() > 0);
        svarname_ << info.svars[j];
        if ( this->traits_.io_type == IOBase<IOType>::GIO_POSIX ) {
          sprintf(cfname_, spformat.str().c_str(), info.odir.c_str(),
                  svarname_.str().c_str(), info.index, myrank);
          fname_.assign(cfname_);
          nb = write_posix(fname_, info, *u[j]); // only writes one variable per file
        }
        else { // GIO_COLL, collective write
          sprintf(cfname_, scformat.str().c_str(), info.odir.c_str(),
                  svarname_.str().c_str(), info.index);
          fname_.assign(cfname_);
          ostate[0] = u[j];
          nb = write_coll(fname_, info, ostate);
        }
        nd = sz_header(info,this->traits_) + u[j]->size()*sizeof(GFTYPE);
        assert(nb == nd && "Incorrect number of bytes written");
      }
    }
    else {                      // multiple components per file
      assert(this->traits_.io_type == IOBase<IOType>::GIO_COLL && "Invalid io_type");
      svarname_.str(""); svarname_.clear();
      assert(filename.length() > 0);
      svarname_ << filename;
      sprintf(cfname_, scformat.str().c_str(), info.odir.c_str(),
              svarname_.str().c_str(), info.index);
      nb = write_coll(fname_, info, u);
      nd = sz_header(info,this->traits_) + u.size()*u[0]->size()*sizeof(GFTYPE);
      assert(nb == nd && "Incorrect number of bytes written");
    }

} // end, write_state_impl


//**********************************************************************************
//**********************************************************************************
// METHOD : read_state_impl
// DESC   : Read files for entire state, based on StateInfo. 
//          
// ARGS   : filepref: filename prefix, if info.multivar > 0; else
//                    info.svars are used to form filenames
//          info    : StateInfo structure
//          u       : state object. Method will attempt to fill all components of u.
//          bstate  : if == TRUE, read state; else read just stateinfo. Default is TRUE.
// RETURNS: none
//**********************************************************************************
template<typename IOType>
void GIO<IOType>::read_state_impl(std::string filepref, StateInfo &info, State  &u, bool bstate)
{
  GString              serr ="read_state_impl: ";
  GINT                 myrank = GComm::WorldRank(comm_);
  GINT                 ivers, nc, nr, nt;
  GElemType            igtype;
  GString              sgtype;
  std::stringstream    format;  // POSIX name format
  std::stringstream    cformat; // collective name format
  State                istate(1);
  Traits               ttraits;

  if ( !bInit_ ) init();

  nc = this->grid_->gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM; 

  if ( bstate && info.sttype == 0 ) { // Check if reading grid components:
    assert(u.size() == nc && "Incorrect number of grid components");
  }

  resize(this->traits_.wfile);
  
  format .str(""); format .clear();
  cformat.str(""); cformat.clear();
  if ( info.sttype ==0 ) { // is a physical state
    format    << "%s/%s.%0" << this->traits_.wtime << "d.%0" << this->traits_.wtask << "d.out";
  }
  else {                   // is a grid file
    format    << "%s/%s.%0" << this->traits_.wtask << "d.out";
  }
  cformat    << "%s/%s.%0" << this->traits_.wtime << "d.out";

  if ( !this->traits_.multivar ) { // one state comp per file

    assert(info.svars.size() >= 1);
    nt = bstate ? u.size() : 1;
    for ( GSIZET j=0; j<nt; j++ ) { // Retrieve all state/grid components
      svarname_.str(""); svarname_.clear();
      assert(info.svars[j].length() > 0);
      svarname_ << info.svars[j];
      sprintf(cfname_, format.str().c_str(), info.idir.c_str(),
              svarname_.str().c_str(), info.index, myrank);
      fname_.assign(cfname_);
      if ( this->traits_.io_type == IOBase<IOType>::GIO_POSIX ) { // POSIX
        nr = read_posix(fname_, info, *u[j], bstate);
      }
      else {  // collective
        istate[0] = u[j];
        nr = read_coll (fname_, info, istate, bstate);
      }
    }

  }
  else {                      // multiple state components in file

    assert(this->traits_.io_type == IOBase<IOType>::GIO_COLL && "Invalid io_type");
    svarname_.str(""); svarname_.clear();
    assert(filepref.length() > 0);
    svarname_ << filepref;
    sprintf(cfname_, format.str().c_str(), info.odir.c_str(),
            svarname_.str().c_str(), info.index);
    read_coll(fname_, info, u, bstate);

  }

} // end, read_state_impl method


//**********************************************************************************
//**********************************************************************************
// METHOD : read_state_info_impl
// DESC   : Read files for StateInfo.
//          
// ARGS   : filepref: filename prefix 
//          info    : StateInfo structure, returned
// RETURNS: none
//**********************************************************************************
template<typename IOType>
void GIO<IOType>::read_state_info_impl(std::string filename, StateInfo &info)
{
  GString              serr ="read_state_info_impl: ";
  GSIZET               nh;
  Traits               ttraits;


  nh  = read_header(filename, info, ttraits);

  assert( nh == sz_header(info,this->traits_) );

} // end, read_state_info_impl method


//**********************************************************************************
//**********************************************************************************
// METHOD : write_posix
// DESC   : Do simple (lowes-level) GIO POSIX output of specified field
// ARGS   : 
//          filename: filename, fully resolved
//          info    : StateInfo structure
//          u       : field to output
// RETURNS: number bytes written
//**********************************************************************************
template<typename IOType>
GSIZET GIO<IOType>::write_posix(GString filename, StateInfo &info, const GTVector<Value> &u) 
{

    GString  serr ="write_posix: ";
    FILE     *fp;
    GSIZET    i, j, nb, nd;

   
    fp = fopen(filename.c_str(),"wb");
    if ( fp == NULL ) {
      cout << serr << "Error opening file: " << filename<< endl;
      exit(1);
    }

    if ( this->traits_.ivers > 0 && info.porder.size(1) <= info.nelems ) {
      cout << serr << " porder of insufficient size for version: " << this->traits_.ivers
                   << ". Error writing file " << filename << endl;
      exit(1);
    }

    // Write header: dim, numelems, poly_order, etc:
    nb = write_header_posix(fp, info, this->traits_);
 
    assert(nb == sz_header(info,this->traits_));

    // Write field data:
    nb += fwrite(u.data(), sizeof(Value), u.size(), fp);
    fclose(fp);

    return nb;

} // end, write_posix


//**********************************************************************************
//**********************************************************************************
// METHOD : read_posix
// DESC   : Do simple GIO POSIX read
// ARGS   : 
//          filename: filename, fully resolved
//          info    : StateInfo structure, filled here
//          u       : field to output; resized if required
//          bstate  : if == TRUE, read state; else read just stateinfo. Default is TRUE.
// RETURNS: none
//**********************************************************************************
template<typename IOType>
GSIZET GIO<IOType>::read_posix(GString filename, StateInfo &info, GTVector<Value> &u, bool bstate)
{

    GString serr = "read_posix: ";
    FILE     *fp;
    GSIZET    nb, nd, nh, nt;
    Traits    ttraits;
    
    nh  = read_header(filename, info, ttraits);

    assert(ttraits.ivers == this->traits_.ivers
                                         && "Incompatible file version number");
    assert(ttraits.dim   == GDIM         && "File dimension incompatible with GDIM");
    assert(info   .gtype == this->grid_->gtype() 
                                           && "File grid type incompatible with grid");
    if ( !bstate ) return nh;

    fp = fopen(filename.c_str(),"rb");
    if ( fp == NULL ) {
      cout << serr << "Error opening file: " << filename << endl;
      exit(1);
    }

    // Seek to first byte after header:
    fseek(fp, nh, SEEK_SET); 

    // Compute field data size from header data:
    nd = 0;
    if ( ttraits.ivers == 0 ) { // expansion order is constant
      nt = 1; 
      for ( GSIZET j=0; j<GDIM; j++ ) nt *= (info.porder(0,j) + 1);
      nd += nt * info.nelems;
    }
    else {                     // expansion order varies
      for ( GSIZET i=0; i<info.porder.size(1); i++ ) {
        nt = 1; 
        for ( GSIZET j=0; j<GDIM; j++ ) nt *= (info.porder(i,j) + 1);
        nd += nt;
      }
    }

    u.resize(nd);
    nb = fread(u.data(), sizeof(Value), nd, fp);

    fclose(fp);

    if ( nb != nd ) {
      cout << serr << "Incorrect amount of data read from file: " << filename << endl;
      exit(1);
    }

   return nb;

} // end, read_posix


//**********************************************************************************
//**********************************************************************************
// METHOD : write_header_posix
// DESC   : Write GIO POSIX file header. Note file pointer position 
//          may be changed on exit.
// ARGS   : 
//          fp       : file pointer (must be opened). Not closed here.
//          info     : StateInfo structure, filled with what header provides
//          traits   : object's traits
// RETURNS: no. header bytes written
//**********************************************************************************
template<typename IOType>
GSIZET GIO<IOType>::write_header_posix(FILE *fp, StateInfo &info, Traits &traits)
{
    GString serr ="write_header: ";
    GSIZET nb, nd, nh, numr;
  
    assert(fp != NULLPTR && "error opening file");

    nb = 0;
  
    fseek(fp, 0, SEEK_SET);
    // Write header: dim, numelems, poly_order:
    nh=fwrite(&traits.ivers     , sizeof(GINT)  ,    1, fp); // GIO version number
      nb += nh*sizeof(GINT);
    nh=fwrite(&traits.dim       , sizeof(GINT)  ,    1, fp); // dimension
      nb += nh*sizeof(GINT);
    nh=fwrite(&info.nelems      , sizeof(GSIZET),    1, fp); // num elements
      nb += nh*sizeof(GSIZET);
    nd = info.porder.size(1) * info.porder.size(2);
    nh=fwrite(info.porder.data().data()
                               , sizeof(GINT)  ,   nd, fp); // exp order
      nb += nh*sizeof(GINT);
    nh=fwrite(&info.gtype         , sizeof(GINT)  ,    1, fp); // grid type
      nb += nh*sizeof(GINT);
    nh=fwrite(&info.cycle         , sizeof(GSIZET),    1, fp); // time cycle stamp
      nb += nh*sizeof(GSIZET);
    nh=fwrite(&info.time          , sizeof(GFTYPE),    1, fp); // time stamp
      nb += nh*sizeof(GFTYPE);
  
    return nb;

} // end, write_header_posix


#if defined(_G_USE_MPI)
//**********************************************************************************
//**********************************************************************************
// METHOD : write_header_coll
// DESC   : Write GIO MPI file header. Note file pointer position 
//          may be changed on exit.
// ARGS   : 
//          fp       : file pointer (must be opened). Not closed here.
//          info     : StateInfo structure, filled with what header provides
//          traits   : object's traits
// RETURNS: no. header bytes written
//**********************************************************************************
template<typename IOType>
GSIZET GIO<IOType>::write_header_coll(MPI_File fp, StateInfo &info, Traits &traits)
{
    GString serr ="write_header: ";
    GINT   nh;
    GSIZET nb, numr;
  

    nb = 0;
  
    MPI_Status     status;
    MPI_File_seek(fp, 0, MPI::SEEK_SET); // set to 0-displacement
    
    MPI_File_write(fp, &traits.ivers   , 1   , T2GCDatatype  <GINT>(), &status); 
        MPI_Get_count(&status, MPI_BYTE, &nh);  nb += nh;
    MPI_File_write(fp, &traits.dim     , 1   , T2GCDatatype  <GINT>(), &status); 
        MPI_Get_count(&status, MPI_BYTE, &nh);  nb += nh;
    MPI_File_write(fp, &info.nelems    , 1   , T2GCDatatype<GSIZET>(), &status); 
        MPI_Get_count(&status, MPI_BYTE, &nh); nb += nh;
    numr = traits.ivers == 0 ? 1 : info.nelems;
    info.porder.resize(numr,traits.dim);
    numr = info.porder.size(1)*info.porder.size(2);
    MPI_File_write(fp, info.porder.data().data()
                                               , numr, T2GCDatatype   <GINT>(), &status); 
        MPI_Get_count(&status, MPI_BYTE, &nh); nb += nh;
    MPI_File_write(fp, &info.gtype       , 1   , T2GCDatatype  <GINT>(), &status);
        MPI_Get_count(&status, MPI_BYTE, &nh);  nb += nh;
    MPI_File_write(fp, &info.cycle       , 1   , T2GCDatatype<GSIZET>(), &status);
        MPI_Get_count(&status, MPI_BYTE, &nh); nb += nh;
    MPI_File_write(fp, &info.time        , 1   , T2GCDatatype<GFTYPE>(), &status); 
        MPI_Get_count(&status, MPI_BYTE, &nh); nb += nh;

    return nb;

} // end, write_header_coll
#endif


//**********************************************************************************
//**********************************************************************************
// METHOD : read_header
// DESC   : Read GIO file header
// ARGS   : 
//          filename : file name (fully resolved)
//          info     : StateInfo structure, filled with what header provides
//          traits   : this object's traits
// RETURNS: no. header bytes read
//**********************************************************************************
template<typename IOType>
GSIZET GIO<IOType>::read_header(GString filename, StateInfo &info, Traits &traits)
{

    GString serr ="read_header: ";
    GSIZET nb, nd, nh, numr;
  
    nb = 0;
//  if ( traits.io_type == IOBase<IOType>::GIO_POSIX ) {
      // Read field data:
      FILE *fp;
      fp = fopen(filename.c_str(),"rb");
      assert(fp != NULLPTR && "gio.cpp: error opening file");
    
      // Read header: 
      nh = fread(&traits.ivers     , sizeof(GINT)  ,    1, fp); nb += nh*sizeof(GINT);
      nh = fread(&traits.dim       , sizeof(GINT)  ,    1, fp); nb += nh*sizeof(GINT);
      nh = fread(&info.nelems      , sizeof(GSIZET),    1, fp); nb += nh*sizeof(GSIZET);
      numr = traits.ivers == 0 ? 1 : info.nelems;
      info.porder.resize(numr,traits.dim);
      numr = info.porder.size(1)*info.porder.size(2);
      nh = fread(info.porder.data().data()
                                   , sizeof(GINT)  , numr, fp); nb += nh*sizeof(GINT);
      nh = fread(&info.gtype       , sizeof(GINT)  ,    1, fp); nb += nh*sizeof(GINT);
      nh = fread(&info.cycle       , sizeof(GSIZET),    1, fp); nb += nh*sizeof(GSIZET);
      nh = fread(&info.time        , sizeof(GFTYPE),    1, fp); nb += nh*sizeof(GFTYPE);
    
      fclose(fp);
  
//  }
/*
    else {
#if defined(_G_USE_MPI)
      MPI_File        fh;
      MPI_Status      status;
      MPI_File_open(comm_, filename.c_str(), MPI::MODE_RDWR, MPI::INFO_NULL, &fh);
//    MPI_File_seek(fh, 0, MPI::SEEK_SET); // set to 0-displacement
      
      // Read header: 
      nh = MPI_File_read(fh, &traits.ivers     , 1   , T2GCDatatype<GINT>  (), &status); nb += nh*sizeof(GINT);
      nh = MPI_File_read(fh, &traits.dim       , 1   , T2GCDatatype<GINT>  (), &status); nb += nh*sizeof(GINT);
      nh = MPI_File_read(fh, &info.nelems      , 1   , T2GCDatatype<GSIZET>(), &status); nb += nh*sizeof(GSIZET);
      numr = traits.ivers == 0 ? 1 : info.nelems;
      info.porder.resize(numr,traits.dim);
      numr = info.porder.size(1)*info.porder.size(2);
      nh = MPI_File_read(fh, info.porder.data().data()
                                               , numr, T2GCDatatype  <GINT>(), &status); nb += nh*sizeof(GINT);
      nh = MPI_File_read(fh, &info.gtype       , 1   , T2GCDatatype<GINT>  (), &status); nb += nh*sizeof(GINT);
      nh = MPI_File_read(fh, &info.cycle       , 1   , T2GCDatatype<GSIZET>(), &status); nb += nh*sizeof(GSIZET);
      nh = MPI_File_read(fh, &info.time        , 1   , T2GCDatatype<GFTYPE>(), &status); nb += nh*sizeof(GFTYPE);

      MPI_File_close(&fh); 
#else
# error "MPI not defined; cannot read GIO_COLL"
#endif
*/
//  } 

    // Check number read vs expected value:
    nd = sz_header(info, traits);
    if ( nb != nd ) {
      cout << serr << "Incorrect amount of data read from file: " << filename << endl;
      exit(1);
    }

    return nb;

} // end, read_header


//**********************************************************************************
//**********************************************************************************
// METHOD : sz_header
// DESC   : Get header size in bytes
// ARGS   : 
//          info     : StateInfo structure, filled with what header provides
//          traits   : object's traits
// RETURNS: no. header bytes 
//**********************************************************************************
template<typename IOType>
GSIZET GIO<IOType>::sz_header(const StateInfo &info, const Traits &traits)
{

    GString serr ="sz_header: ";
    GSIZET nd, numr;
  
    // Get no. bytes in header (should agree with read_posx, write_posix):
    numr = traits.ivers == 0 ? 1 : info.nelems;
    numr *= traits.dim;
    nd = (numr+3)*sizeof(GINT) + 2*sizeof(GSIZET) + sizeof(GFTYPE);

    return nd;

} // end, sz_header


//**********************************************************************************
//**********************************************************************************
// METHOD : resize
// DESC   : Resize global data
// ARGS   : n  : new number bytes
// RETURNS: none.
//**********************************************************************************
template<typename IOType>
void GIO<IOType>::resize(GINT n)
{

  if ( n > nfname_ ) {
    if ( cfname_ != NULLPTR ) delete [] cfname_;
    cfname_ = new char [n];
    fname_.resize(n);
    nfname_ = n;
  }

} // end, resize


//**********************************************************************************
//**********************************************************************************
// METHOD : write_coll
// DESC   : Collective write of state components
// ARGS   : filename: filename
//          info    : StateInfo structure
//          u       : state
// RETURNS: number bytes written
//**********************************************************************************
template<typename IOType>
GSIZET GIO<IOType>::write_coll(GString filename, StateInfo &info, const State &u)
{
#if !defined(_G_USE_MPI)
  #error "Illegal entry into GIO<IOType>::write_coll: MPI not defined"
#endif

#if defined(_G_USE_MPI)

    GString        serr = "write_coll: ";
    GINT           myrank = GComm::WorldRank(comm_);
    GINT           iret, nc, nh;
    GSIZET         nb;
    GSIZET         ntot;

    MPI_Offset     disp;
    MPI_File       fh;
    MPI_Status     status;

    // Required number of coord vectors:
    nc = this->grid_->gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;


    iret = MPI_File_open(comm_, filename.c_str(), MPI::MODE_CREATE|MPI::MODE_WRONLY, MPI::INFO_NULL, &fh);
    assert(iret == MPI_SUCCESS && "MPI_File_open failure");
    nbheader_ = sz_header(info, this->traits_);

    nh = write_header_coll(fh, info, this->traits_);
    assert(nh == nbheader_ && "Expected header size not written");

    ntot = nh;

    // Cycle over all fields, and write:
    if ( !this->traits_.multivar ) { // print each comp to sep. file
      for ( auto j=0; j<u.size(); j++ ) {
        assert(u[j]->size() > 0 && "Invalid state component");
        disp = nbheader_ + j*nbgdof_ + state_disp_;
        iret = MPI_File_set_view(fh, disp, T2GCDatatype<GFTYPE>(), mpi_state_type_, "native", MPI::INFO_NULL);
        iret = MPI_File_write(fh, u[j]->data(), 1, mpi_state_type_, &status);
        MPI_Get_count(&status, MPI_BYTE, &nh);  
        ntot += nh;
      }
    }
    else {                       // print each comp to same file
      for ( auto j=0; j<u.size(); j++ ) {
        disp = nbheader_ + state_disp_;
        iret = MPI_File_set_view(fh, disp, T2GCDatatype<GFTYPE>(), mpi_state_type_, "native", MPI::INFO_NULL);
        iret = MPI_File_write(fh, u[j]->data(), 1, mpi_state_type_, &status);
        MPI_Get_count(&status, MPI_BYTE, &nh);  
        ntot += nh;
      }
    } 

    MPI_File_close(&fh);
#endif

    return ntot;

} // end, write_coll


//**********************************************************************************
//**********************************************************************************
// METHOD : read_coll
// DESC   : Collective read of state components
// ARGS   : filename: filename
//          info    : StateInfo structure
//          u       : state
//          bstate  : if == TRUE, read state; else read just stateinfo. Default is TRUE.
// RETURNS: none
//**********************************************************************************
template<typename IOType>
GSIZET GIO<IOType>::read_coll(GString filename, StateInfo &info, State &u, bool bstate)
{
#if defined(_G_USE_MPI)

    GString        serr = "read_coll: ";
    GINT           myrank = GComm::WorldRank(comm_);
    GINT           iret, nc;
    GSIZET         nb, nh, ntot;
    Traits         ttraits;
    StateInfo      tinfo;

    MPI_Offset     disp;
    MPI_File       fh;
    MPI_Status     status;

    // Required number of coord vectors:
    nc = this->grid_->gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;


    MPI_File_open(comm_, filename.c_str(), MPI::MODE_RDONLY, MPI::INFO_NULL, &fh);
    nbheader_ = sz_header(info, this->traits_);

    // Read header and do some checks:
    nh = read_header(filename, tinfo, ttraits);
    assert(nh == nbheader_ && "Expected header size not read");
    assert(ttraits.ivers   == this->traits_.ivers && "Incompatible GIO version");   
    assert(ttraits.dim     == this->traits_.dim   && "Incompatible problem dimension");
    assert(tinfo.porder      == info.porder    && "Incompatible porder");
    assert(tinfo.gtype       == info.gtype     && "Incompatible grid type");

    ntot = nh;

    if ( !bstate ) return nh;

    // Cycle over all fields, and read:
    //   Note: any variable order element-by-element
    //         should be handled in ::init method via 
    //         state_disp_ & mpi_state_type_. Currently,
    //         variable order is not fully supported on read:
    if ( !this->traits_.multivar ) { // print each comp to sep. file
      for ( auto j=0; j<u.size(); j++ ) {
        assert(u[j]->size() > 0 && "Invalid state component");
        disp = nbheader_ + j*nbgdof_ + state_disp_;
        iret = MPI_File_set_view(fh, disp, T2GCDatatype<GFTYPE>(), mpi_state_type_, NULL, MPI::INFO_NULL);
        iret = MPI_File_read(fh, u[j]->data(), 1, mpi_state_type_, &status);
        ntot += iret == MPI_SUCCESS ? u[j]->size() : 0;
      }
    }
    else {                       // read each comp from same file
      for ( auto j=0; j<u.size(); j++ ) {
        disp = nbheader_ + state_disp_;
        iret = MPI_File_set_view(fh, disp, T2GCDatatype<GFTYPE>(), mpi_state_type_, NULL, MPI::INFO_NULL);
        iret = MPI_File_read(fh, u[j]->data(), 1, mpi_state_type_, &status);
        ntot += iret == MPI_SUCCESS ? u[j]->size() : 0;
      } 
    } 

    MPI_File_close(&fh);
#endif

    return ntot;

} // end, read_coll

