function x = roundn(x,n)
%ROUNDN rounds x to n didgets.
%
% USAGE:
%    x = roundn(x,n);
%
% TO 121108

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

n=round(n);

x = round(x*10.^n)/10.^n;