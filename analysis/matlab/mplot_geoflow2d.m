function h = mplot_geoglow2d(svar, tindex)
%
% Does a mesh plot of 2D GeoFLOW data
%
%  Usage:
%    h = mplot_geoflow2d('u1',10)
%
%  Input:
%    s1var   : prefix for field file. Required
%    tindex  : time index for output. Required
%
%  Output:
%    h       : plot handle
%
if nargin < 2
  error('must specify svar and tindex');
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


  % Cycle over elems, and plot 'patches':
  icurr = 1;
  for n = 1:nelems
    xx = x{1}(icurr:icurr+lelem-1);
    yy = x{2}(icurr:icurr+lelem-1);
    uu = u   (icurr:icurr+lelem-1);
    xx = reshape(xx, NN(1), NN(2));
    yy = reshape(yy, NN(1), NN(2));
    uu = reshape(uu, NN(1), NN(2));
    surf( xx, yy, uu, 'edgecolor', 'none' )
%   surf( xx, yy, uu)
    title(sprintf('%s t=%f', svar, time));
    hold on
    icurr = icurr + lelem ; 
    
  end % end, elem loop
  

end % end, task loop

