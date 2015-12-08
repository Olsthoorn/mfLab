function Idx = outside(xm,a,b)
%OUTSIDE Get indices if xm outside a given range [a b];
%
% Example:
%   Idx = outside(xm,a,b)
%
% See also:  hit between fallsIn above below after beyond outside inMesh
%
%   TO 1204020

Idx = find(find(xm<min(a,b) | xm>max(a,b)));
      