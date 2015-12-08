function matlabdate=excel2datenum(exceldate)
%EXCEL2DATENUM converts excel datenum to matlab datenum
%
% USAGE:
%    matlabdate=excel2datenum(exceldate)
%
% Copyright TO 2009-1012

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

matlabdate = exceldate + datenum(1899,12,30);
