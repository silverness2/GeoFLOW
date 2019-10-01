//==================================================================================
// Module       : gtstat.ipp
// Date         : 10/1/19 (DLR)
// Description  : Encapsulates the methods and data associated with
//                GeoFLOW statistics methods.
// Copyright    : Copyright 2019, Colorado State University. All rights reserved.
// Derived From : none.
//==================================================================================



//**********************************************************************************
//**********************************************************************************
// METHOD : Constructor method
// DESC   : Instantitate with no. bins
// ARGS   : GSIZET n: no. bins
// RETURNS: none
//**********************************************************************************
template<class T>
GTStat<T>::GTStat(GSIZET n, GC_COMM comm):
nbins_         (n),
nkeep_         (0),
comm_       (comm)
{
  assert(nbins_ > 0 && "Invalid bin count");

  myrank_  = GComm::WorldRank(comm_);

  lpdf  .resize(nbins_);
  ikeep_.resizem(u.size();
}


//**********************************************************************************
//**********************************************************************************
// METHOD : Copy constructor method
// DESC   : Override default copy constructor
// ARGS   : 
// RETURNS: 
//**********************************************************************************
template<class T>
GTStat<T>::GTStat(const GTStat<T> &obj):
nbins_      (obj.get_nbins())
{
  lpdf  .resize(nbins_);
  ikeep_.resizem(u.size();
}


//**********************************************************************************
//**********************************************************************************
// METHOD : dopdf1d (1)
// DESC   : Do 1d pdf of input scalar field, and output to specified file.
// ARGS   : 
//          u      : scalar field on which to operate
//          ifixdr : if TRUE, then use specified fmin, fmax to set dynamical range.
//                   If FALSE, then determine the dynamical range dynamically, and
//                   provide these to caller
//          fmin/
//          fmax   : dynamical range for pdf
//          dolog  : take log of |u| when creating bins?
//          pdf    : final pdf
// RETURNS: none.
//**********************************************************************************
template<class T> 
void dopdf1d(GTVector<T> u, GBOOL ifixdr, T &fmin, T &fmax, GDOOL dolog, GTVector<T> &pdf)
{
  T      del, test, tmax, tmin;
  T      fbin, fmin1, fmax1;
  T      gavg, sig, sumr, xnorm;
  T      tiny;
  GSIZET ibin, lkeep;

  gpdf_.resizem(nbins_);
  lpdf_ = 0.0;
  pdf   = 0.0;

  tiny = fabs(100.0*std::numeric_limits<T>::epsilon());

  // Compute dynamic range, if not using specified range:
  if ( !ifixdr ) {
    fmin1 = u.min();
    fmax1 = u.min();
    GComm::Allreduce(fmin1, fmin, 1, T2GCDatatype<T>() , GC_OP_MIN, comm_);
    GComm::Allreduce(fmax1, fmax, 1, T2GCDatatype<T>() , GC_OP_MAX, comm_);
  }
  
  if ( dolog ) {
    fmin = log10(fabs(fmin)+tiny);
    fmax = log10(fabs(fmax)+tiny);
  }

  tmin = fmin;
  tmax = fmax;
  if ( dolog .GT. 0 ) {
    tmin = pow(10.0, fmin);
    tmax = pow(10.0, fmax);
  }

  // Find indirection indices that meet dyn. range criterion:
  lkeep = 0;
  #pragma omp for
  for ( auto j=0; j<u.size(); j++ ) {
    if ( u[j] >= tmin && u[i] <= tmax ) {
      ikeep_[lkeep] = j;
      lkeep++;
    }
  }
  GComm::Allreduce(lkeep, nkeep_, 1, T2GCDatatype<GSISET>() , GC_OP_SUM, comm_);
  assert(nkeep_ > 0  && "No samples are within dynamic range");

  xnorm = 1.0 / std::static_cast<T>(nkeep);

  // Compute average:
  sumr = 0.0;
  #pragma omp parallel for default(shared) private(j) reduction(+:sumr)
  for ( auto j=0; j<nkeep_; j++ ) {
    sumr += u[ikeep_[j]];
  }
  GComm::Allreduce(sumr, gavg, 1, T2GCDatatype<T>() , GC_OP_SUM, comm_);
  gavg *= xnorm;

  // Compute std deviation:
  sumr = 0.0;
  #pragma omp parallel for default(shared) private(j) reduction(+:sumr)
  for ( auto j=0; j<nkeep_; j++ ) {
    sumr += pow(u[ikeep_[j]]-gavg,2);
  }
  GComm::Allreduce(sumr, gavg, 1, T2GCDatatype<T>() , GC_OP_SUM, comm_);
  sig *= xnorm;

  // Compute local PDF:
  del = fabs(fmax - fmin) / nbins;
  if ( dolog ) {
    #pragma omp parallel for  private(ibin,test)
    for ( auto j=0; j<nkeep_; j++ ) {
      test = log10(fabs(u(ikeep_[i]))+tiny);
      ibin = std::static_cast<GSIZET> ( ( test - fmin )/del+1 );
      ibin = MIN(MAX(ibin,1),nbins_);
      #pragma omp atomic
      lpdf_[ibin] += 1.0;
    }
  }
  else {
    #pragma omp parallel for  private(ibin,test)
    for ( auto j=0; j<nkeep_; j++ ) {
      test = fabs(u(ikeep_[i]));
      ibin = std::static_cast<GSIZET> ( ( test - fmin )/del+1 );
      ibin = MIN(MAX(ibin,1),nbins_);
      #pragma omp atomic
      lpdf_[ibin] += 1.0;
    }
  }
  
 // Do sanity check:
  fbin = 0.0;
  #pragma omp parallel for  default(shared) reduction(+:fbin)
  for ( auto j=0; j<nbins_; j++ ) {
    fbin += pdf_[j]; 
    #pragma omp atomic
    fbin += lpdf_(ibin);
  }
  assert( fbin == nkeep_ && "Inconsistent binning");

  // Compute global reduction between MPI tasks to find final (global) pdf:
  GComm::Allreduce(lpdf_, pdf, nbins_, T2GCDatatype<T>() , GC_OP_SUM, comm_);


} // end, dopdf1d (1)


//**********************************************************************************
//**********************************************************************************
// METHOD : dopdf1d (2)
// DESC   : Do 1d pdf of input scalar field, and return pdf. Additional member
//          is allocated in this call (gpdf_).
// ARGS   : 
//          u      : scalar field on which to operate
//          ifixdr : if TRUE, then use specified fmin, fmax to set dynamical range.
//                   If FALSE, then determine the dynamical range dynamically, and
//                   provide these to caller
//          fmin/
//          fmax   : dynamical range for pdf
//          dolog  : take log of |u| when creating bins?
//          fname  : filename to which to output pdf
// RETURNS: none.
//**********************************************************************************
template<class T> 
void dopdf1d(GTVector<T> u, GBOOL ifixdr, T &fmin, T &fmax, GDOOL dolog, const GString &fname)
{
  T                 ofmax, ofmin;
  std::ofstream     ios;
  std::stringstream header;

  gpdf_.resizem(npdf_);

  dopdf1d(u, ifixdr, fmin, fmax, dolog, gpdf_);
  if ( myrank.eq.0 )  {

    if ( dolog ) {
      ofmin = pow(10.0, fmin);
      ofmax = pow(10.0, fmax);
    }
    WRITE(shead,'(A9,E16.8,A1,E16.8,A7,E16.8,A6,E16.8,A7,I7,A7,I1,A8,F12.0)') '# range=[', &
        fmin, ',' , fmax, ']; avg=',gavg,'; sig=',sig,'; nbin=', nbins, '; blog=', dolog,   &
        '; nkeep=',gkeep
    header << std::scientific << setprecision(8)
    header << "# range=[" 
           << ofmin  << ","
           << ofmax  << "]; avg="
           << gavg   << "; sig="
           << sig    << "; nbins="
           << nbins_ << "; blog ="
           << dolog  << "; nkeep="
           << nkeep_ << std::endl;

     ios.open(fname,std::ios_base::trunc);    
     ios << header << std:endl 
     ios << std::scientific << setprecision(15)
     for ( auto j=0; j<nbins_-1; j++ ) {
       ios << gpdf_[j] << " " 
     }
     ios << gpdf_[nbins_-1] << std::endl;
    
     ios.close();
  }


} // end, dopdf1d (2)


