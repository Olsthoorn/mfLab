function Idx = beyond(xm,a)
%BEYOND get indices of cells xGr or xm beyond a where xGr}xm is assumed ascending
%
% Example:
%    Idx = beyond(xm,a)
%
% See also:  hit between fallsIn above below after beyond outside inMesh
%
%   TO 1204020

   Idx= find(xm>a);
   
   