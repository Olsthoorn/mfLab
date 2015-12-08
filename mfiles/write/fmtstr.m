function str=fmtstr(fmt,ntimes)
%FMTSTR writes ntimes a format (i.e. %2d or so) to a string for use in sprintf
%
% Example:
%    str=fmtstr(fmt,ntimes);
%
% TO 090101

fmt=fmt(ones(ntimes,1),:)';

str=['',fmt(:)','\r\n'];   % eg ' %i%i%i %i %i\n'
