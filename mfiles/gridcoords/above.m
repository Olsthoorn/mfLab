function Idz = above(zGr,a)
%ABOVE gets indices of zGr or zm before point a, where zGr is ascending
%
% Example:
%   Idz = above(xGr,a)
%
% See also:  hit between fallsIn above below after beyond outside inMesh
%
%   TO 1204020

Idz = before(zGr,a);
