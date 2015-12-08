function bool = isaxis(var)
%ISAXIS verify that var is an axis
%
% USAGE:
%    logical = isaxis(variable)
%
%    Notice that axis must be child of a figure.
%
% TO120531 130331 130402

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

bool = ishghandle(var,'axes');
