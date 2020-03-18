function h = boxplot2d(svar, tindex, blog, dtype, ftype, varargin)
%
% Does a line plot 2D GeoFLOW Posix data in x-direction, 
% at fixed specified y-index
%
%  Usage:
%    h = boxplot2d('u1',10, [0 1])
%
%  Input:
%    svar    : prefix for field file. Required
%    tindex  : time index for output. Required
%    blog    : take log of data?
%    bwire   : do wire plot?
%    dtype   : data file type: 'POSIX', 'COLLC'. Default is 'COLL'
%    ftype   : floating point size (4 or 8). Default is 8.
%    varargin: to pass to quadmesh: e.g. to plot
%              wire mesh only and set to single color,
%              set bwire=1, and call:
%
%              boxplot2d('u1', 1, 0, 1, 'POSIX', 8, 'edgecolor','b')
%                 or
%              boxplot2d('u1', 1, 0, 0, 'COLL' , 8, 'colorbarlims',[-2 2])
%
%  Output:
%    h       : plot handle
%
if nargin < 2
  error('must specify svar and tindex');
end 

if nargin < 3
  blog = 0;
  dtype = 'COLL';;
  ftype = 8;
end 
if nargin < 4
  dtype = 'COLL';
  ftype = 8;
end 
if nargin < 5
  ftype = 8;
end 

if dtype ~= 'POSIX' && dtype ~= 'COLL'
  error(['Invalid dtype: ' dtype]);
end

sz = size(plottype)
if ~isempty(plottype) && sz(1)*sz(2) ~= 2 
  error('incorrect plottype specification');
end

lwidth = 2;
szfont = 16;

vartmp = varargin;
bcolorbarlims = 0;
n = length(varargin);
j = 1;
while j <= n
  if strcmpi(vartmp{j},'colorbarlims') == 1
    if n < j+1
      error(sprintf('colorbarlims variable requires array [a b] of limits'));
    end
    colorbarlims = vartmp{j+1}
    bcolorbarlims = 1
    varargin = { vartmp{1:j-1}, vartmp{j+2:end} }
    j = j + 1
  end
  j = j + 1
end

scoord = {'xgrid','ygrid' 'zgrid'};

ntasks = 1;
if dtype == 'POSIX'
  d = dir('xgrid.*');
  ntasks = length(d);
  if ntasks<= 0 
    error('Grid data missing or incomplete');
  end
end

% Find global size:
if dtype == 'POSIX'
  tsize = zeros(ntasks,1); % total size per task
  for itask = 0:ntasks-1
    fname = sprintf('%s.%06d.%05d.out', svar, tindex, itask);
    [u dim nelems porder gtype icycle time] = rgeoflow(fname, ftype, 'ieee-le', 1);
    NN = double(porder + 1); 
  % lelem = prod(NN(1:dim));  % data length per element
    tsize (itask+1) = NN(1)*nelems;
  end
  nglobal = sum(tsize); % global no. nodes
else
  fname = sprintf('%s.%06d.out', svar, tindex);
  [u dim nelems porder gtype icycle time] = rgeoflow(fname, ftype, 'ieee-le', 1);
  NN = double(porder + 1); 
end


igstart = 1;
iegsave = 1;
ntot    = 0;
for itask = 0:ntasks-1

  % Find start index in global data for this 
  % task's data:
  if itask > 0
    igstart = sum(tsize(1:itask)) + 1;
  end

  % Read node coords:
  for j=1:2
    if dtype == 'POSIX'
      fname = sprintf('%s.%05d.out', scoord{j}, itask);
    elseif dtype == 'COLL'
      fname = sprintf('%s.out', scoord{j});
    end
    [x{j} dim nelems porder gtype icycle time mvar] = rgeoflow(fname, ftype, 'ieee-le');
  end
  if ( dim ~= 2 )
    error('Grid must be 2D');
  end 

  if dtype == 'POSIX'
    fname = sprintf('%s.%06d.%05d.out', svar, tindex, itask);
  elseif dtype == 'COLL'
    fname = sprintf('%s.%06d.out', svar, tindex);
  end
  [u dim nelems porder gtype icycle time mvar] = rgeoflow(fname, ftype, 'ieee-le');

 
  NN = double(porder + 1); 
  lelem = prod(NN(1:dim));  % data length per element

  if jindex < 0 || jindex >= NN(2) 
    error('Invalid jindex');
  end

  % Cycle over elems, and gather 22 data:
  icurr = 1;
  for n = 1:nelems
    xx  = x{1}(icurr:icurr+lelem-1);
    yy  = x{2}(icurr:icurr+lelem-1);
    uu  = u   (icurr:icurr+lelem-1);
    pdorder = double(porder);

    % Find indices defining quads between
    % node points:
    imat = zeros(int64(prod(pdorder)/nverts), nverts);
    n = 1;
    for k = 1:pdorder(2) % cycle over each quad
      for j = 1:pdorder(1)
        bb = [ j+(k-1)*NN(1) j+1+(k-1)*NN(1) j+1+k*NN(1) j+k*NN(1) ];
        imat(n,:) = bb;
        n = n + 1;
      end
    end
    cu = uu(imat);
    cm = mean(cu,2);
    cf = (cm - umin) / (umax - umin + eps);
    cv = (uu - umin) / (umax - umin + eps);
    if bwire == 0
      if blog > 0
        puu = log10(uu);
      else
        puu = uu;
      end
      h = quadmesh(imat,xx,yy,zz,puu,'FaceColor','interp');
      colorbar('vertical');
    else
      h = quadmesh(imat,xx,yy,zz,varargin{:});
    end
    hold on;
%   set(p, 'FaceColor', 'blue', 'EdgeColor', 'none');
    if bwire == 0
      title(sprintf('%s t=%f', svar, time));
    end


    icurr = icurr + lelem; 
    
  end % end, elem loop

end % end, task loop


end
