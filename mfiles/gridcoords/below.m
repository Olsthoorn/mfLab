function Idz = below(zGr,z0)
%BELOW  gets indices of zGr or zm before point a, where zGr/zm are ascending
%
% Example:
%    Idx = below(zGr,a)
%
% See also:  hit between fallsIn above below after beyond outside inMesh
%
%   TO 1204020

Idz = after(zGr,z0);   