function Idx = fallsIn(xGr,a,b)
%FALLSIN  gets indices indices of xGr or xm between find(xGr>min(a,b) & xGr<max(a,b)).
%
% Example:
%   Idx = fallsIn(xGr,a,b);
%   Idx = between(xGr,[a b]);
%
%   Indices outside the range b a if b<a;
%
% See also:  hit between fallsIn above below after beyond outside inMesh
%
%   TO 1204020  130421

if nargin<2, error('%s: minimum 2 inputs',mfilename); end

if nargin==2, b=a(end); a=a(1); end

i1 = hit(xGr,a);
i2 = hit(xGr,b);
Idx= min(i1,i2):max(i1,i2);

   