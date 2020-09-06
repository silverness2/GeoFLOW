function [c h] = boxcont2d(svar, tindex, blog, xlev, levtype, nozero, dtype, isz )
%
%   function h = boxcont2d(svar, tindex, blog, xlev, levtype, nozero, dtype, isz)
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
%    xlev    : number of levels, if levtype = 'num'; contour
%              increment if levtype = 'delta'. Default is 5 if
%              not specified, and levtype='num' is the default.
%    levtype : if 'num', then lev is the number of levels to use.
%              If levtype=='delta', then lev is the difference between
%    nozero  : if >0, don't include zero-contour; else do. Default is 1.
%              contour levels. Default is 'num'.
%    dtype   : data file type: 'POSIX', 'COLL'. Default is 'COLL'
%    isz    : floating point size (4 or 8). Default is 8.
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
  xlev  = 5;
  levtype = 'num';
  nozero = 0;
  dtype = 'COLL';;
  isz = 8;
end 
if nargin < 4
  xlev  = 5;
  levtype = 'num';
  nozero = 0;
  dtype = 'COLL';
  isz = 8;
end 
if nargin < 5
  levtype = 'num';
  nozero = 0;
  dtype = 'COLL';
  isz = 8;
end 
if nargin < 6
  nozero = 0;
  dtype = 'COLL';
  isz = 8;
end 
if nargin < 7
  dtype = 'COLL';
  isz = 8;
end 
if nargin < 8
  isz = 8;
end 

dtype
if ~strcmp(dtype,'POSIX') & ~strcmp(dtype,'COLL')
  error(['Invalid dtype: ' dtype]);
end


lwidth = 2;
szfont = 16;

bcolorbarlims = 0;

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


  dx = diff(x{1}(:));
  I  = find(abs(dx) > 0.0 );
  ngridx = ( max(x{1}(:)) - min(x{1}(:)) ) / min(abs(dx(I)))
  ngridy = ngridx;
  
  [X,Y] = ndgrid(linspace(min(x{1}(:)),max(x{1}(:)),ngridx),linspace(min(x{2}(:)),max(x{2}(:)),ngridy));

  Z = griddata(x{1}(:),x{2}(:),u(:),X,Y,'cubic');

Z(200:300)

  zmin = min(min(Z));
  zmax = max(max(Z));
  if strcmpi(levtype,'num')
    nlev = int32(xlev);
    dc = (zmax - zmin)/(double(nlev)+1);
    dcvec = [0:double(nlev)].*dc + zmin;
  elseif strcmpi(levtype,'delta')
    nlev = int32((zmax - zmin)/xlev);
    dc = (zmax - zmin)/(double(nlev)+1);
    dcvec = [0:double(nlev)].*dc + zmin;
  else
    error(['Invalid levtype: ' levtype]);
  end

  % Remove zero-contour:
  if nozero > 0
    II = find(abs(dcvec) > 1.0e-6);
    levels = dcvec(II);
  else
    levels = dcvec;
  end

  [c h] = contour(X, Y, Z, levels); 
  axis equal

  hold on;

end % end, task loop

title(sprintf('%s t=%f', svar, time));
xlabel('x (m)');
ylabel('z (m)');

end
