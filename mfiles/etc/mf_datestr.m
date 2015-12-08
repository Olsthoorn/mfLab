function s = mf_datestr(d,dateFormat,format)
%MF_DATESTR sprints a data using matlab datestr if >31/12/1799 and format if otherwise.
%
% Usage:
%    s = mf_datestr(d,dateFormat)
%
%    if d>datenum(1799,12,31), same as datestr(d,dateformat)
%    else s = sprintf(dateFormat,d)
%
% Examples:
%   s = mf_datestr(now,'yyyy-mmm-dd');
%   s = mf_datestr( 5.2,'%g');
%
% See also: datestr
%
% TO 130408


if d>datenum(1799,12,31)
    s = datestr(d,dateFormat);
else
    if nargin<=2
        format=' %.4g';
    end
    s = sprintf('%.4g %s',d,format);
end
    