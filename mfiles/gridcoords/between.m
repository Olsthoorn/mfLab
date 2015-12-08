function Idx = between(xGr,a,b)
%BETWEEN  gets indices indices of xGr or xm between find(xGr>min(a,b) & xGr<max(a,b)).
%
% Example:
%   Idx = between(xm,a,b);
%   Idx = between(xm,[a b]);
%   
%   Indices outside the range b a if b<a;
%
% See also:  hit between fallsIn above below after beyond outside inMesh
%
%   TO 1204020

if nargin<2, error('%s: minimum 2 inputs',mfilename); end

if nargin==2, b=a(end); a=a(1); end

Idx = find(xGr>min(a,b) & xGr<max(a,b));

   