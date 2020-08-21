function [c h] = boxcont2d(svar, tindex, blog, nlev, dtype, isz, varargin)
%
%   function h = boxcont2d(svar, tindex, blog, nlev, dtype, isz, varargin)
%
% Plots contour plots for 2D GeoFLOW Posix or collective data. 
%
%  Usage:
%    [c h] = boxcont2d('u1',10, 0, 8)
%
%  Input:
%    svar    : prefix for field file. Required
%    tindex  : time index for output. Required
%    blog    : take log of data? Default is 0
%    nlev    : number of levels. Default is 5
%    dtype   : data file type: 'POSIX', 'COLL'. Default is 'COLL'
%    isz    : floating point size (4 or 8). Default is 8.
%    varargin: 
%
%  Output:
%    c       : contour matrix
%    h       : plot handle
%
if nargin < 2
  error('must specify svar and tindex');
end 

if nargin < 3
  blog = 0;
  nlev = 5;
  dtype = 'COLL';;
  isz = 8;
end 
if nargin < 4
  nlev = 5;
  dtype = 'COLL';
  isz = 8;
end 
if nargin < 5
  dtype = 'COLL';
  isz = 8;
end 
if nargin < 6
  isz = 8;
end 

dtype
if ~strcmp(dtype,'POSIX') & ~strcmp(dtype,'COLL')
  error(['Invalid dtype: ' dtype]);
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

[umin, umax] = gminmax_gio(svar, tindex, dtype, isz, 'ieee-le');

ntasks = 1;
if strcmp(dtype,'POSIX')
  d = dir('xgrid.*');
  ntasks = length(d);
  if ntasks<= 0 
    error('Grid data missing or incomplete');
  end
end

nverts = 4;

figure();
for itask = 0:ntasks-1

  % Read node coords:
  for j=1:2
    if strcmp(dtype,'POSIX')
      fname = sprintf('%s.%06d.%05d.out', scoord{j}, 0, itask);
    elseif dtype == 'COLL'
      fname = sprintf('%s.%06d.out', scoord{j}, 0);
    end
    [x{j} dim nelems porder gtype icycle time mvar] = rgeoflow(fname, isz, 'ieee-le');
  end
  if ( dim ~= 2 )
    error('Grid must be 2D');
  end 
dim
nelems

  if strcmp(dtype,'POSIX')
    fname = sprintf('%s.%06d.%05d.out', svar, tindex, itask);
  elseif dtype == 'COLL'
    fname = sprintf('%s.%06d.out', svar, tindex);
  end
  [u dim nelems porder gtype icycle time mvar] = rgeoflow(fname, isz, 'ieee-le');

 
  NN = double(porder + 1); 
  if blog == 1 
    u = log10(abs(u));
  end

  [X,Y] = ndgrid(linspace(min(x{1}(:)),max(x{1}(:)),100),linspace(min(x{2}(:)),max(x{2}(:)),100));

  Z = griddata(x{1}(:),x{2}(:),u(:),X,Y,'cubic');
X(100:200)
Y(200:300)
Z(200:300)
  [c h] = contour(X, Y, Z, nlev); 
  axis equal

end % end, task loop


end
