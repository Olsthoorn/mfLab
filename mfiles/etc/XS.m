function A=XS(A)
%XS convert 3D array to cross section along X-axis?
%
% USAGE:
%    XS A=XS(A)
%
% See also YS
%
% TO 090101 120410
%

% Copyright 2009 2012 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if nargin<1,  error('Use A = XS(A) or A=YS(A) for permutations to face ZX and ZY planes respectively'); end

A=permute(A,[3 2 1]);

end
