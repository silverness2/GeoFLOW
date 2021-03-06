function [xg  ug h] = boxplot1d(svar, tindex, blog, symb, jindex)
%
% Does a line plot 2D GeoFLOW Posix data in x-direction, 
% at fixed specified y-index
%
%  Usage:
%    [xg ug h] = boxplot1d('u1',10, [0 1])
%
%  Input:
%    s1var   : prefix for field file. Required
%    tindex  : time index for output. Required
%    blog    : array with 2 elems; a 1 in one element
%              implies log on correspoing axis. 
%              I.e., [1 0] means log-lin, [0 1] lin-log, 
%              [1 1], log-log. If empty, no plot is done.
%    dtype   : data file type: 'POSIX', 'COLL'. Default is 'COLL'
%    symb    : plot symbol. Default is 'k-' (blaxk solid)
%    jindex  : j-index. Default is 0.
%
%  Output:
%    xg      : global x coord
%    ug      : global field at each x
%    h       : plot handle
%
if nargin < 2
  error('must specify svar and tindex');
end 
if nargin < 3
 blog = [0 0];
 dtype = 'COLL';
 symb  = 'k-';
 jindex   = 1;
end 
if nargin < 4
 dtype = 'COLL';
 symb     = 'k-';
 jindex   = 1;
end 
if nargin < 5
 symb     = 'k-';
 jindex   = 1;
end 
if nargin < 6
 jindex   = 1;
end 
sz = size(blog)
if ~isempty(blog) && sz(1)*sz(2) ~= 2 
  error('incorrect glob specification');
end

if ~strcmp(dtype,'POSIX') & ~strcmp(dtype,'COLL')
  error(['Invalid dtype: ' dtype]);
end

doplot = 1
if isempty(blog)
  doplot = 0
end

lwidth = 2;
szfont = 16;

ntasks = 2;
scoord = {'xgrid','ygrid' 'zgrid'};

ntasks = 1;
if strcmp(dtype,'POSIX')
  d = dir('xgrid.*');
  ntasks = length(d);
  if ntasks<= 0
    error('Grid data missing or incomplete');
  end
end

% Find global size:
tsize = zeros(ntasks,1); % total size per task
for itask = 0:ntasks-1
  if strcmp(dtype,'POSIX')
    fname = sprintf('%s.%06d.%05d.out', svar, tindex, itask);
  elseif dtype == 'COLL'
    fname = sprintf('%s.%06d.out', svar, tindex);
  end

  [u dim nelems porder gtype icycle time mvar] = rgeoflow(fname, 8, 'ieee-le', 1);
  NN = double(porder + 1); 
% lelem = prod(NN(1:dim));  % data length per element
  tsize (itask+1) = NN(1)*nelems;
end
nglobal = sum(tsize); % global no. nodes


% Allocate global data:
xg = zeros(nglobal,1);
ug = zeros(nglobal,1);

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
    if strcmp(dtype,'POSIX')
      fname = sprintf('%s.%06d.%05d.out', scoord{j}, 0, itask);
    elseif dtype == 'COLL'
      fname = sprintf('%s.%06d.out', scoord{j}, 0);
    end
    [x{j} dim nelems porder gtype icycle time mvar] = rgeoflow(fname, 8, 'ieee-le');
  end
  if ( dim ~= 2 )
    error('Grid must be 2D');
  end 

  if strcmp(dtype,'POSIX')
    fname = sprintf('%s.%06d.%05d.out', svar, tindex, itask);
  elseif dtype == 'COLL'
    fname = sprintf('%s.%06d.out', svar, tindex);
  end

  [u dim nelems porder gtype icycle time mvar] = rgeoflow(fname, 8, 'ieee-le');

 
  NN = double(porder + 1); 
  lelem = prod(NN(1:dim));  % data length per element

  if jindex < 0 || jindex >= NN(2) 
    error('Invalid jindex');
  end

  % Cycle over elems, and gather 1d data:
  icurr = 1;
  for n = 1:nelems
    xx  = x{1}(icurr:icurr+lelem-1);
%   yy  = x{2}(icurr:icurr+lelem-1);
    uu  = u   (icurr:icurr+lelem-1);
    ib  = jindex*NN(1) + 1; % index local to elem 
    ie  = ib + NN(1) - 1;   % index local to elem
     
    ibg = igstart + (n-1)*NN(1);     % beg index in global array
    ieg = ibg + NN(1) - 1;           % end index in global array

    if ieg > nglobal
      error('...................inconsistent data');
    end
%   disp(sprintf('ibg=%d; ieg=%d',ibg,ieg))

    xg(ibg:ieg) = xx(ib:ie);
    ug(ibg:ieg) = uu(ib:ie);
    iegsave = ieg;

    icurr = icurr + lelem; 
    
    ntot = ntot + NN(1);
  end % end, elem loop

end % end, task loop


if ntot ~= size(xg)
ntot
size(xg)
  error('total processed does not match');
end

% Plot 1d profile:
[xg I] = sort(xg,1);
ug    = ug(I);
h = -1;
if doplot > 0
  figure
  if ( blog(1) == 0 && blog(2) == 0 )
    h = plot(xg, ug, symb);
  elseif ( blog(1) == 1 && blog(2) == 0 )
    h = semilogx(abs(xg), ug, symb);
  elseif ( blog(1) == 0 && blog(2) == 1 )
    h = semilogy(xg, abs(ug), symb);
  elseif ( blog(1) == 1 && blog(2) == 1 )
    h = loglog(abs(xg), abs(ug), symb);
  else
    error('invalid blog');
  end
end

if doplot > 0

  %title(sprintf('%s index=%d  t=%f', svar, tindex, time));
  title(sprintf('index=%d  t=%f', tindex, time));
  xlabel('x','FontSize',szfont,'FontWeight','Bold');
  ylabel(strtok(svar),'FontSize',szfont,'FontWeight','Bold');
  set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',lwidth,'FontSize',szfont,'FontWeight','Bold');
  
end
