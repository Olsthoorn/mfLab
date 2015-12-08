function a = ismynan(X,dim)
%ISMYNAN returns true of any X==NaN
%
% USAGE:
%    a = ismynan(X,dim)
% 
% same along dim if dim is given 
%
%    usefull in cellfun ot make sure that only a single output is produced per
%    cell
%
% See also cellfun

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if nargin<2
    a = any(isnan(X));
else
    a = any(isnan(X,dim));
end