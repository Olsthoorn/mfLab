function Idx = after(xGr,a)
%AFTER Indices of cells after a where xGr or xm is assumed to ascend
%
% Example:
%    Idx = after(xGr,a)
%
% See also:  hit between fallsIn above below after beyond outside inMesh
%
%    TO 1204020

Idx=find(xGr>a);
      
   