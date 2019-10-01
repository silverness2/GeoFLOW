function h = boxplot1d(svar, tindex, jindex)
%
% Does a line plot 2D GeoFLOW data in x-direction, at fixed specified y-index
%
%  Usage:
%    h = boxplot1d('u1',10)
%
%  Input:
%    s1var   : prefix for field file. Required
%    tindex  : time index for output. Required
%    jindex  : j-index. Default is 0.
%
%  Output:
%    h       : plot handle
%
if nargin < 2
  error('must specify svar and tindex');
end 
if nargin < 3
  jindex = 0;
end 


ntasks = 2;
scoord = {'xgrid','ygrid' 'zgrid'};

d = dir('xgrid.*');
ntasks = length(d);
if ntasks<= 0 
  error('Grid data missing or incomplete');
end

for itask = 0:ntasks-1

  % Read node coords:
  for j=1:2
    fname = sprintf('%s.%05d.out', scoord{j}, itask)
    [x{j} dim nelems porder gtype icycle time] = rgeoflow(fname, 8, 'ieee-le');
  end
  if ( dim ~= 2 )
    error('Grid must be 2D');
  end 

  fname = sprintf('%s.%06d.%05d.out', svar, tindex, itask);
  [u dim nelems porder gtype icycle time] = rgeoflow(fname, 8, 'ieee-le');
  if ( itask == 0 )
    figure;
  end

 
  NN = double(porder + 1); 
  lelem = prod(NN(1:dim))  % data length per element

  if jindex < 0 || jindex >= NN(2) 
    error('Invalid jindex');
  end

  xn = zeros(NN(1)*nelems,1);
  un = zeros(NN(1)*nelems,1);

  % Cycle over elems, and gather 1d data:
  icurr = 1;
  for n = 1:nelems
    xx = x{1}(icurr:icurr+lelem-1);
%   yy = x{2}(icurr:icurr+lelem-1);
    uu = u   (icurr:icurr+lelem-1);
    ib = jindex*NN(1) + 1 
    ie = ib + NN(1) 
     
    ibn = (n-1)*NN(1) + 1 
    ien = ibn + NN(1) 
    xn(ibn:ien) = xx(ib:ie);
    un(ibn:ien) = uu(ib:ie);

    icurr = icurr + lelem ; 
    
  end % end, elem loop
  
  % Plot 1d profile:
  [xn I] = sort(xn);
   un    = un(I);
  plot(xn, un, 'k-');

end % end, task loop

