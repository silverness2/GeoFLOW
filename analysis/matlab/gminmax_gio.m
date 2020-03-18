function [umin umax] = gminmax_gio(svar, tindex, dtype, isz, sformat)
%
% Finds global max, min of specified variable
%
%  Usage:
%    [dmin, dmax]  = gminmax_gio('u1',10)
%
%  Input:
%    svar    : prefix for field file. Required
%    tindex  : time index for output. Required
%    dtype   : write-type: 'POSIX' or 'COLL'. Default is 'COLL'
%    isz     : data size (in bytes: either 4 or 8, e.g.). Default is isz=8.
%    sformat : data format of file: 'ieee-be' or 'ieee-le' for big-endian or little
%              endian if isz=4, or 'ieee-be.l64', 'ieee-le.l64' if isz=8. Default
%              is 'ieee-le'.

%
%  Output:
%    umax   
%    umin   : max, min vals

%

if nargin < 5
  error('must specify svar, tindex, isz, sformat');
end 
if nargout < 2
  error('Too few output parameters');
end 

if ~strcmp(dtype,'POSIX') & ~strcmp(dtype,'COLL')
  error(['Invalid dtype: ' dtype]);
end


ntasks = 1;
if strcmp(dtype,'POSIX')
  d = dir(sprintf('%s.%06d.*',svar, tindex));
  ntasks = length(d);
  if ntasks<= 0 
    error('Grid data missing or incomplete');
  end
end

umax   = -realmax;
umin   =  realmax;

for itask = 0:ntasks-1

  if strcmp(dtype,'POSIX')
    fname = sprintf('%s.%06d.%05d.out', svar, tindex, itask)
  elseif strcmp(dtype,'COLL')
    fname = sprintf('%s.%06d.out', svar, tindex)
  end
  [u dim nelems porder gtype icycle time mvar] = rgeoflow(fname, isz, sformat);
 
  umin = min(umin, min(u));
  umax = max(umax, max(u));

end % end, task loop


end
