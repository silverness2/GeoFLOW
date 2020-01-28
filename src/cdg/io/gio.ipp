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
IOBase<IOType>(),
bInit_                     (FALSE),
comm_                       (comm),
default_state_name_pref_ ("state"),
default_grid_name_pref_   ("grid"),
cfname_                  (NULLPTR),
nfname_                        (0)
{ 
  traits_  = &traits;
  grid_    = &grid;
#if !defined(_G_USE_MPI)
  assert(trants_->io_type == GIO_COLL && "Collective IO only allowed if MPI is used");
#endif

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
  MPI::Type::free( &mpi_state_type_);
#endif
} // end of destructor method


//**********************************************************************************
//**********************************************************************************
// METHOD     : write_state_impl
// DESCRIPTION: Write state to file
//
// ARGUMENTS  : filename: file to write to
//              info    : state info structure
//              u       : state to write
//               
// RETURNS    : none.
//**********************************************************************************
template<typename IOType>
void GIO<IOType>::write_state_impl(std::string filename, StateInfo &info, State &u)
{

  switch ( traits_->io_type ) {

    case GIO_POSIX:
      write_state_posix(info, u);
      break;
    case GIO_COLL:
      write_state_coll(info, u);
      break;
    default :
      assert(FALSE);
  }

} // end of method write_state_impl


//**********************************************************************************
//**********************************************************************************
// METHOD     : read_state_impl
// DESCRIPTION: Read state from file
//
// ARGUMENTS  : filename: file to read from to
//              info    : state info structure
//              u       : storage for state
//               
// RETURNS    : none.
//**********************************************************************************
template<typename IOType>
void GIO<IOType>::read_state_impl(std::string filename, StateInfo &info, State &u)
{

  switch ( traits_->io_type ) {

    case GIO_POSIX:
      read_state_posix(info, u);
      break;
    case GIO_COLL:
      read_state_coll(info, u);
      break;
    default :
      assert(FALSE);
  }


} // end of method read_state_impl


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

  if ( traits_->io_type != GIO_COLL ) {
    bInit_ = TRUE;
    return; // nothing to do
  }

  GSIZET           ndof  =grid_->ndof();
  GSIZET           nelems=grid_->nelems();
  GTVector<GSIZET> extent;

  extent.resize(GComm::WorldSize(comm_));
 
#if defined(_G_USE_MPI)


  GComm::Allgather(&ndof, 1, T2GCDatatype<GSIZET>, extent.data(), 1, T2GCDatatype<GSIZET>, comm_);

  GINT myrank = GComm::WorldRank(comm_);
  if ( myrank == 0 ) {
    state_disp_ = 0;
    state_extent_ = extent[myrank];
  }
  else {
    state_disp_   = extent.sum(0,myrank)*sizeof(GFTYPE);
    state_extent_ = extent[myrank]*sizeof(GFTYPE);
  }
  nbgdof_ = extent.sum();

  MPI::Aint lowerbnd=0;
  MPI::Type::create_resized(T2GCDatatype<GFTYPE>, lowerbnd, state_extent_, &mpi_state_type_);
  MPI::Type::commit(&mpi_state_type_);

#endif


} // end of method init



//**********************************************************************************
//**********************************************************************************
// METHOD : write_state_posix
// DESC   : Do simple GIO POSIX output of state. Each state member is
//          written to its own file. Each file is tagged s.t.:
//             svar.TTTTTT.PPPPP.out,
//          where TTTTTT is the time tag, and PPPPP is the task number. 
//          svar is the name of each state field, provided in info structre. 
//          If there aren't enough svar specified, then defaults are used.
//          Note: if info.sttype > 0, then we assume we are printing a 
//          grid, and the format is the same as above but without the time 
//          index tag.
// ARGS   : info: StateInfo structure
//          u     : state
// RETURNS: none
//**********************************************************************************
template<typename IOType>
void GIO<IOType>::write_state_posix(StateInfo &info, const State &u)
{
    GString        serr = "write_state_posix: ";
    GINT           myrank = GComm::WorldRank(comm_);
    GSIZET         nb, nc, nd;
    std::vector<GString> 
                   gpref = {"x","y","z"};
    GTVector<GTVector<GFTYPE>>
                  *xnodes = &grid_->xNodes();
    GElemList     *elems  = &grid_->elems();
    std::stringstream 
                   format;

    // Required number of coord vectors:
    nc = grid_->gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;
    if ( info.sttype !=0 ) { // is a grid file 
      assert(nc ==  xnodes->size()
          && "Incompatible grid or coord  dimensions");
    }
    
    traits.dim  = GDIM;
    info.nelems = grid_->nelems(); // local nelems for this task

    info.gtype    = grid_->gtype();


    // Build format string:
    resize(traits.wfile);
    if ( info.sttype ==0 ) { // is a physical state
      format    << "%s/%s.%0" << traits.wtime << "d.%0" << traits.wtask << "d.out";
    }
    else {                   // is a grid file
      format    << "%s/%s.%0" << traits.wtask << "d.out";
    }

    // Set porder vector depending on version:
    info.porder.resize(traits.ivers == 0 ? 1 : traits.nelems, GDIM);
    for ( auto i=0; i<info.porder.size(1); i++ )  { // for each element
      for ( auto j=0; j<info.porder.size(2); j++ ) info.porder(i,j) = (*elems)[i]->order(j);
    }

    // Cycle over all fields, and write:
    svarname_.str("");
    svarname_.clear();
    for ( auto j=0; j<u.size(); j++ ) {
      assert(u[j].size() > 0 && "Invalid state component");
      if ( info.svars.size() < u.size() ||  into.svars[j].length() <= 0 ) {
        if ( info.sttype ==0 ) { // is a physical state
          svarname_ << default_state_name_pref_ << j;
        }
        else {
          svarname_ << gpref[j] << default_grid_name_pref_;
        }
      }
      sprintf(cfname_, format.str().c_str(), info.odir.c_str(),
              svarname_.str().c_str(), traits.index, myrank);
      fname_.assign(cfname_);
      write_posix<Value>(fname_, traits, *u[j]);
    }

} // end, write_state_posix


//**********************************************************************************
//**********************************************************************************
// METHOD : read_state_posix
// DESC   : Read restart files for entire state, based on StateInfo. Since
//          the state components reside in different files, each state element
//          is read using as file name  prefixes from the StateInfo.svars labels. 
//          The fully-resolved path is created by prepending the info.idir, and 
//          appending the task index.         
//          
// ARGS   : 
//          info  : StateInfo structure
//          u     : state object. Method will attempt to fill all components of u.
// RETURNS: none
//**********************************************************************************
template<typename IOType>
void GIO<IOType>::read_state_posix(StateInfo &info, State  &u)
{
  GString              serr ="read_state_posix: ";
  GINT                 myrank = GComm::WorldRank(comm_);
  GINT                 ivers, nc, nr, nt;
  GElemType            igtype;
  GSIZET               itindex;
  GString              sgtype;
  std::stringstream    format;
  Traits               objtraits;
  StateInfo            info;

  nc = grid_->gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM; 

  if ( info.sttype == 0 ) { // Check if reading grid components:
    assert(u.size() == nc && "Incorrect number of grid components");
  }

  resize(traits.wfile);
  format    << "%s/%s.%0" << traits_->wtask << "d.out";
  for ( GSIZET j=0; j<u.size(); j++ ) { // Retrieve all state/grid components
    if ( into.svars[j].length() <= 0 ) {
      cout << serr << "empty filename for component " << j << endl; 
      exit(1);
    }
    sprintf(cfname_, format.str().c_str(), info.idir.c_str(), svars[j].c_str(),  myrank);
    fname_.assign(cfname_);
    nr = read_posix<GFTYPE>(info.svars[j], info, *u[j]);
  }

} // end, read_state_posix method


//**********************************************************************************
//**********************************************************************************
// METHOD : write_posix
// DESC   : Do simple (lowes-level) GIO POSIX output of specified field
// ARGS   : 
//          filename: filename, fully resolved
//          info    : StateInfo structure
//          u       : field to output
// RETURNS: total number bytes written
//**********************************************************************************
template<typename IOType>
GSIZET GIO<IOType>::write_posix(GString filename, StateInfo &info, const GTVector<T> &u) 
{

    GString  serr ="write_posix: ";
    FILE     *fp;
    GSIZET    i, j, nb, nd;

   
    fp = fopen(filename.c_str(),"wb");
    if ( fp == NULL ) {
      cout << serr << "Error opening file: " << filename << endl;
      exit(1);
    }

    if ( traits_->ivers > 0 && info.porder.size(1) <= info.nelems ) {
      cout << serr << " porder of insufficient size for version: " << traits_->ivers
                   << ". Error writing file " << filename << endl;
      exit(1);
    }

    // Write header: dim, numelems, poly_order, etc:
    nb = write_header<FILE*>(fp, info);
 
    assert(nb == sz_header(info,*traits_));

    // Write field data:
    nb += fwrite(u.data(), sizeof(T), u.size(), fp);
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
// RETURNS: none
//**********************************************************************************
template<typename IOType>
GSIZET GIO<IOType>::read_posix(GString filename, StateInfo &info, GTVector<T> &u)
{

    GString serr = "read_posix: ";
    FILE     *fp;
    GSIZET    nb, nd, nh, nt;
    Traits    objtraits;
    
    nh  = read_header(filename, info, objtraits);

    assert(objtraits.ivers == traits_->ivers
                                           && "Incompatible file version number");
    assert(objtraits.dim   == GDIM         && "File dimension incompatible with GDIM");
    assert(info     .gtype == grid_->gtype() 
                                           && "File grid type incompatible with grid");

    fp = fopen(filename.c_str(),"rb");
    if ( fp == NULL ) {
      cout << serr << "Error opening file: " << fname << endl;
      exit(1);
    }

    // Seek to first byte after header:
    fseek(fp, nh, SEEK_SET); 

    // Compute field data size from header data:
    nd = 0;
    if ( objtraits.ivers == 0 ) { // expansion order is constant
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
    nb = fread(u.data(), sizeof(T), nd, fp);

    fclose(fp);

    if ( nb != nd ) {
      cout << serr << "Incorrect amount of data read from file: " << fname << endl;
      exit(1);
    }

   return nb;

} // end, read_posix


//**********************************************************************************
//**********************************************************************************
// METHOD : write_header
// DESC   : Write GIO file header. Note file pointer position 
//          may be changed on exit.
// ARGS   : 
//          fp       : file pointer (must be opened). Not closed here.
//          info     : StateInfo structure, filled with what header provides
// RETURNS: no. header bytes written
//**********************************************************************************
template<typename IOType>
template<typename GFPTR>
GSIZET GIO<IOType>::write_header(GFPTR fp, StateInfo &info)
{
    GString serr ="write_header: ";
    GSIZET nb, nd, nh;
  
    assert(fp != NULLPTR && "error opening file");

    nb = 0;
  
    if ( traits_->io_type == GIO_POSIX ) {
    fseek(fp, 0, SEEK_SET);
    // Write header: dim, numelems, poly_order:
    nh=fwrite(&traits_->ivers     , sizeof(GINT)  ,    1, fp); // GIO version number
      nb += nh*sizeof(GINT);
    nh=fwrite(&traits_->dim       , sizeof(GINT)  ,    1, fp); // dimension
      nb += nh*sizeof(GINT);
    nh=fwrite(&info.nelems        , sizeof(GSIZET),    1, fp); // num elements
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
  }
  else {

#if defined(_G_USE_MPI)
      MPI::Status     status
      MPI::File::seek(fh, 0, MPI::SEEK_SET); // set to 0-displacement
      
      nh = MPI::File::write(fh, &traits.ivers     , 1   , T2GCDatatype  <GINT>(), &status); nb += nh*sizeof(GINT);
      nh = MPI::File::write(fh, &traits.dim       , 1   , T2GCDatatype  <GINT>(), &status); nb += nh*sizeof(GINT);
      nh = MPI::File::write(fh, &info.nelems      , 1   , T2GCDatatype<GSIZET>(), &status); nb += nh*sizeof(GSIZET);
      numr = traits.ivers == 0 ? 1 : info.nelems;
      info.porder.resize(numr,traits.dim);
      numr = info.porder.size(1)*info.porder.size(2);
      nh = MPI::File::write(fh, info.porder.data().data()
                                               , numr, T2GCDatatype   <GINT>(), &status); nb += nh*sizeof(GINT);
      nh = MPI::File::write(fh, &info.gtype       , 1   , T2GCDatatype  <GINT>(), &status); nb += nh*sizeof(GINT);
      nh = MPI::File::write(fh, &info.cycle       , 1   , T2GCDatatype<GSIZET>(), &status); nb += nh*sizeof(GSIZET);
      nh = MPI::File::write(fh, &info.time        , 1   , T2GCDatatype<GFTYPE>(), &status); nb += nh*sizeof(GFTYPE);

#endif
  }
  
    return nb;

} // end, write_header


//**********************************************************************************
//**********************************************************************************
// METHOD : read_header
// DESC   : Read GIO file header
// ARGS   : 
//          filename : file name (fully resolved)
//          info     : StateInfo structure, filled with what header provides
//          traits   : object's traits
// RETURNS: no. header bytes read
//**********************************************************************************
template<typename IOType>
GSIZET GIO<IOType>::read_header(GString filename, StateInfo &info, Traits &traits)
{

    GString serr ="read_header: ";
    GSIZET nb, nd, nh, numr;
  
    nb = 0;
    if ( traits_->io_type == GIO_POSIX ) {
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
  
    }
    else {
#if defined(_G_USE_MPI)
      MPI::File       fh;
      MPI::Status     status
      MPI::File::open(comm_, fname_.c_str(), MPI::MODE_RDWR, MPI::INFO_NULL, &fh);
//    MPI::File::seek(fh, 0, MPI::SEEK_SET); // set to 0-displacement
      
      // Read header: 
      nh = MPI::File::read(fh, &traits.ivers     , 1   , T2GCDatatype<GINT>  (), &status); nb += nh*sizeof(GINT);
      nh = MPI::File::read(fh, &traits.dim       , 1   , T2GCDatatype<GINT>  (), &status); nb += nh*sizeof(GINT);
      nh = MPI::File::read(fh, &info.nelems      , 1   , T2GCDatatype<GSIZET>(), &status); nb += nh*sizeof(GSIZET);
      numr = traits.ivers == 0 ? 1 : info.nelems;
      info.porder.resize(numr,traits.dim);
      numr = info.porder.size(1)*info.porder.size(2);
      nh = MPI::File::read(fh, info.porder.data().data()
                                               , numr, T2GCDatatype  <GINT>(), &status); nb += nh*sizeof(GINT);
      nh = MPI::File::read(fh, &info.gtype       , 1   , T2GCDatatype<GINT>  (), &status); nb += nh*sizeof(GINT);
      nh = MPI::File::read(fh, &info.cycle       , 1   , T2GCDatatype<GSIZET>(), &status); nb += nh*sizeof(GSIZET);
      nh = MPI::File::read(fh, &info.time        , 1   , T2GCDatatype<GFTYPE>(), &status); nb += nh*sizeof(GFTYPE);

      MPI::File::close(fh); 
#endif

    } 

    // Check number read vs expected value:
    nd = (numr+3)*sizeof(GINT) + 2*sizeof(GSIZET) + sizeof(GFTYPE);
    if ( nb != nd ) {
      cout << serr << "Incorrect amount of data read from file: " << fname << endl;
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
GSIZET GIO<IOType>::sz_header(StateInfo &info, Traits &traits)
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
// RETURNS: none
//**********************************************************************************
template<typename IOType>
void GIO<IOType>::write_coll(GString filename, StateInfo &info, const State &u)
{
#if defined(_G_USE_MPI)

    GString        serr = "write_coll: ";
    GINT           myrank = GComm::WorldRank(comm_);
    GSIZET         nb, nh;

    MPI::Offset    disp;
    MPI::File      fh;
    MPI::Status    status

    // Required number of coord vectors:
    nc = grid_->gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;


    MPI::File::open(comm_, fname_.c_str(), MPI::MODE_CREATE|MPI::MODE_WRONLY, MPI::INFO_NULL, &fh);
    nbheader_ = sz_header(info, *traits_);

    nh = write_header<MPI::File>(fh, info);
    assert(nh == nbheader_ && "Expected header size not written");

    // Cycle over all fields, and write:
    if ( !traits_->multivar ) { // print each comp to sep. file
      for ( auto j=0; j<u.size(); j++ ) {
        assert(u[j].size() > 0 && "Invalid state component");
        disp = nbheader_ + j*nbgdof_ + state_disp_;
        MPI::File::set_view(fh, disp, T2GCDatatype<GFTYPE>, mpi_state_type_, NULL, MPI::INFO_NULL);
        MPI::File::write(fh, u[j]->data(), 1, mpi_state_type_, &status);
      }
    }
    else {                       // print each comp to same file
      disp = nbheader_ + state_disp_;
      MPI::File::set_view(fh, disp, T2GCDatatype<GFTYPE>, mpi_state_type_, NULL, MPI::INFO_NULL);
      MPI::File::write(fh, u[j]->data(), 1, mpi_state_type_, &status);
    } 

    MPI::File::close(&fh);
#endif

} // end, write_coll


//**********************************************************************************
//**********************************************************************************
// METHOD : read_coll
// DESC   : Collective read of state components
// ARGS   : filename: filename
//          info    : StateInfo structure
//          u       : state
// RETURNS: none
//**********************************************************************************
template<typename IOType>
void GIO<IOType>::read_coll(GString filename, StateInfo &info, const State &u)
{
#if defined(_G_USE_MPI)

    GString        serr = "read_coll: ";
    GINT           myrank = GComm::WorldRank(comm_);
    GSIZET         nb, nh;

    MPI::Offset    disp;
    MPI::File      fh;
    MPI::Status    status

    // Required number of coord vectors:
    nc = grid_->gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;


    MPI::File::open(comm_, fname_.c_str(), MPI::MODE_RDONLY, MPI::INFO_NULL, &fh);
    nbheader_ = sz_header(info, *traits_);

    nh = read_header<MPI::FILE>(fh, info);
    assert(nh == nbheader_ && "Expected header size not written");


    // Cycle over all fields, and write:
    if ( !traits_->multivar ) { // print each comp to sep. file
      for ( auto j=0; j<u.size(); j++ ) {
        assert(u[j].size() > 0 && "Invalid state component");
        disp = nbheader_ + j*nbgdof_ + state_disp_;
        MPI::File::set_view(fh, disp, T2GCDatatype<GFTYPE>, mpi_state_type_, NULL, MPI::INFO_NULL);
        MPI::File::read(fh, u[j]->data(), 1, mpi_state_type_, &status);
      }
    }
    else {                       // print each comp to same file
      disp = nbheader_ + state_disp_;
      MPI::File::set_view(fh, disp, T2GCDatatype<GFTYPE>, mpi_state_type_, NULL, MPI::INFO_NULL);
      MPI::File::read(fh, u[j]->data(), 1, mpi_state_type_, &status);
    } 

    MPI::File::close(&fh);
#endif

} // end, read_coll

