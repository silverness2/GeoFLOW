//==================================================================================
// Module       : gposixio_observer.ipp
// Date         : 3/18/19 (DLR)
// Description  : Observer object for carrying out simple POSIX-based  
//                binary output.
// Copyright    : Copyright 2019. Colorado State University. All rights reserved
// Derived From : IOBaseType.
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
} // end of constructor (1) method


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

  if ( bInit_ ) return;

  bInit_ = TRUE;

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
    GTVector<GTVector<GFTYPE>>
                  *xnodes = &grid_->xNodes();
    GElemList     *elems  = &grid_->elems();
    std::stringstream 
                   format;

    // Required number of coord vectors:
    nc = grid_->gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;

    // Do some checks:
    assert(svars.size() >=  iu.size()
        && "Insufficient number of state variable names specified");

    
    traits.dim  = GDIM;
    info.nelems = grid_->nelems(); // local nelems for this task

    info.gtype    = grid_->gtype();

    assert(nc ==  xnodes->size()
        && "Incompatible grid or coord  dimensions");

    // Build format string:
    resize(traits.wfile);
    format    << "%s/%s.%0" << traits.wtime << "d.%0" << traits.wtask << "d.out";

    // Set porder vector depending on version:
    info.porder.resize(traits.ivers == 0 ? 1 : traits.nelems, GDIM);
    for ( auto i=0; i<info.porder.size(1); i++ )  { // for each element
      for ( auto j=0; j<info.porder.size(2); j++ ) info.porder(i,j) = (*elems)[i]->order(j);
    }

    // Cycle over all fields, and write:
    svarname_.str("");
    svarname_.clear();
    for ( auto j=0; j<u.size(); j++ ) {
      if ( info.svars.size() < u.size() ||  into.svars[j].length() <= 0 ) {
        svarname_ = default_state_name_pref_ << j;
      }
      sprintf(cfname_, format.str().c_str(), info.odir.c_str(),
              svarname_.str().c_str(), traits.index, myrank);
      assert(u[j].size() > 0 && "Invalid state component index");
      fname_.assign(cfname_);
      write_posix<Value>(fname_, traits, *u[j]);
    }

} // end, write_state_posix


//**********************************************************************************
//**********************************************************************************
// METHOD : write_grid_posix
// DESC   : Do simple GIO POSIX output of grid. Each grid component is
//          taken to be a 'state' compoment, and is written to its own 
//          file
//             svar.PPPPP.out,
//          where PPPPP is the task number. Note: there is no time tag, as
//          grid is assumed not to change.  svar is the name of each 
//          grid component field, provided in info structre. 
//          If there aren't enough svar specified, then defaults are used.
// ARGS   : info : StateInfo structure
// RETURNS: none
//**********************************************************************************
template<typename IOType>
void GIO<IOType>::write_grid_posix(StateInfo &info) 
{
    GString serr ="write_grid_posix: ";
    GINT           myrank = GComm::WorldRank(comm_);
    GSIZET         nb, nc, nd;
    std::vector<GString> 
                   gpref = {"x","y","z"};
    GTVector<GTVector<GFTYPE>>
                  *xnodes = &grid_->xNodes();
    GElemList     *elems  = &grid_->elems();
    std::stringstream format;

    // Required number of coord vectors:
    nc = grid_->gtype() == GE_2DEMBEDDED ? GDIM+1 : GDIM;

    // Do some checks:
    assert(nc ==  xnodes->size()
        && "Incompatible grid or coord  dimensions");

    traits.dim    = GDIM;
    info  .nelems = grid.nelems();
    info  .gtype  = grid_->gtype();


    // Build format string:
    resize(traits.wfile);
    format    << "%s/%s.%0" << traits.wtask << "d.out";
    // Set porder vector depending on version:
    traits.porder.resize(traits.ivers == 0 ? 1: info.nelems, GDIM);
    for ( auto i=0; i<traits.porder.size(1); i++ )  { // for each element
      for ( auto j=0; j<traits.porder.size(2); j++ ) traits.porder(i,j) = (*elems)[i]->order(j);
    }

    // Cycle over all coords, and write:
    sgridname_.str("");
    sgridname_.clear();
    for ( auto j=0; j<xnodes->size(); j++ ) {
      if ( info.svars.size() < u.size() ||  into.svars[j].length() <= 0 ) {
        sgridname_ = gpref[j] << default_grid_name_pref_;
      }
      sprintf(cfname_, format.str().c_str(), info.odir.c_str(), svars[j].c_str(),  myrank);
      fname_.assign(cfname_);
      write_posix<GFTYPE>(fname_, info, (*xnodes)[j]);
    }

#if defined(_G_DEBUG) 
    // Write multiplicity:
    sprintf(cfname_, format.str().c_str(), info.odir.c_str(), "mult",  myrank);
    fname_.assign(cfname_);
    write_posix<GFTYPE>(fname_, info, grid.get_ggfx().get_mult());
 
    GPP(comm_,
                "mult_min="  << grid.get_ggfx().get_mult().min() 
             << " mult_max=" << grid.get_ggfx().get_mult().max() );   
    if ( GComm::WorldRank(comm_) == 0 ) {
      cout << "write: mult=" << grid.get_ggfx().get_mult() << endl;
    }

    // Write glob indices:
    sprintf(cfname_, format.str().c_str(), info.odir.c_str(), "glob_index",  myrank);
    fname_.assign(cfname_);
    write_posix<GNODEID>(fname, info, grid.get_ggfx().glob_index_);
    
#endif

} // end, write_grid_posix


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
// RETURNS: none
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

    // Write header: dim, numelems, poly_order:
    fwrite(&traits_->ivers     , sizeof(GINT)  ,    1, fp); // GIO version number
    fwrite(&traits_->dim       , sizeof(GINT)  ,    1, fp); // dimension
    fwrite(&info.nelems        , sizeof(GSIZET),    1, fp); // num elements
    nd = info.porder.size(1) * info.porder.size(2);
    fwrite(info.porder.data().data()
                               , sizeof(GINT)  ,   nd, fp); // exp order
    fwrite(&info.gtype         , sizeof(GINT)  ,    1, fp); // grid type
    fwrite(&info.cycle         , sizeof(GSIZET),    1, fp); // time cycle stamp
    fwrite(&info.time          , sizeof(GFTYPE),    1, fp); // time stamp

    // Write field data:
    nb = fwrite(u.data(), sizeof(T), u.size(), fp);
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
  
    // Read field data:
    FILE *fp;
    nb = 0;
    fp = fopen(filename.c_str(),"rb");
    assert(fp != NULLPTR && "gio.cpp: error opening file");
  
    // Read header: 
    nh = fread(&traits.ivers       , sizeof(GINT)  ,    1, fp); nb += nh*sizeof(GINT);
    nh = fread(&traits.dim         , sizeof(GINT)  ,    1, fp); nb += nh*sizeof(GINT);
    nh = fread(&info.nelems        , sizeof(GSIZET),    1, fp); nb += nh*sizeof(GSIZET);
    numr = traits.ivers == 0 ? 1 : info.nelems;
    info.porder.resize(numr,traits.dim);
    numr = info.porder.size(1)*info.porder.size(2);
    nh = fread(info.porder.data().data()
                                 , sizeof(GINT)  , numr, fp); nb += nh*sizeof(GINT);
    nh = fread(&info.gtype       , sizeof(GINT)  ,    1, fp); nb += nh*sizeof(GINT);
    nh = fread(&info.cycle       , sizeof(GSIZET),    1, fp); nb += nh*sizeof(GSIZET);
    nh = fread(&info.time        , sizeof(GFTYPE),    1, fp); nb += nh*sizeof(GFTYPE);
  
    fclose(fp);
  
    // Get no. bytest that should have been read:
    nd = (numr+3)*sizeof(GINT) + 2*sizeof(GSIZET) + sizeof(GFTYPE);
  
    if ( nb != nd ) {
      cout << serr << "Incorrect amount of data read from file: " << fname << endl;
      exit(1);
    }

    return nb;

} // end, read_header


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


