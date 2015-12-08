function s = mf_datestr(d,dateFormat,format)
%MF_DATESTR if d>datenum(1799,12,31), same as datestr(d,dateformat)
%
% USAGE:
%    s = mf_datestr(d,dateFormat)
%
% if    d>datenum(1799,12,31), same as datestr(d,dateformat)
% else s = sprintf(dateFormat,d)
%
% Usefull for setting dates in titles.
%
% Example:
%   s = mf_datestr(now,'yyyy-mmm-dd');
%   s = mf_datestr(datenum(1600,36,12));
%   s = mf_datestr(datenum(1600,36,12),'yyy-mmm-dd','%g');
%
% TO 130408

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later


if d>datenum(1799,12,31)
    s = datestr(d,dateFormat);
else
    if nargin<=2
        format=' %.4g';
    end
    s = sprintf('%.4g %s',d,format);
end
    