function exceldate=datenum2excel(matlabdate)
%DATENUM2EXCEL  converts matlab datenum to excel datenum
%
% USAGE:
%    excelDatenum = datenum2excel(matlabDatenum)
%
% TO 110813

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

exceldate = matlabdate - datenum(1899,12,30);
