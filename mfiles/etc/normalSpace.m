function [x,dx]=normalSpace(xMin,xMid,xMax,N,dxRatio,sigma)
%NORMALSPAC generates a nice spacing based on normal distribution
%
% Example:
%   [x,dx]=normalSpace(xMin,xMid,xMax,N,dxRatio,sigma)
%    xMin and xMax determing the extent of the range
%    xMid determines location where dx is dxMin
%    dxRatio = dxMin/dxMax
%    dxMax determines maximum dx at infinite and -infinity
%    sigmal scales x axis as in the normal distribution.
%    notices that xMid does not need to be between xMin and xMax
%    but xMax must be > xMin, dxRation must be >0, N>0 and sigma>0.
%
% See also: linespace logspace logspace2 gridObj gridTutorial
%
% TO 100214 100918


if nargin<6,
    error('Six arguments required: x1,mu,x2,dxMin,dxMax,sigma');
end

if ~(xMax>xMin)
    error('xMax must be larger than xMin');
end

if dxRatio>0
    dxRatio = 1/dxRatio;
end

if ~(sigma>0)
    error('sigma must be > 0');
end

N = round(N);

if~(N>0)
    error('N must be > 0');
end


% xmu 1s center between successive node points
xmu = (xMin:(xMax-xMin)/N:xMax) - xMid; xmu = 0.5*(xmu(1:end-1)+xmu(2:end));

% use normal distribution to compute dx
dx = dxRatio + (1-dxRatio) * (1 - exp(-(xmu/(sigma*sqrt(2))).^2));
dx = dx*(xMax-xMin)/sum(dx);
x  = xMin+[0 cumsum(dx)];
