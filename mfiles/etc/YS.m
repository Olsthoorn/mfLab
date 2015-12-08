function A=YS(A)
%YS convert 3D array to cross section along Y-axis?
%
% USAGE:
%    YS A=XS(A)
%
% see also: XS
%
% Copyright 2009 2012 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

% TO 090101 120410

if nargin<1, error('Use A=YS(A) to permute A onto ZY plane for cross sections along Y'); end

A=permute(A,[3 1 2]);

end
