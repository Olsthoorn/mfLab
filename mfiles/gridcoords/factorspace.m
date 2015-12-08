function y = factorspace(L, fac, n)
%FACTORSPACE get a power series from 0 to L with factor fac and n points
%
%   Example:
%      y = factorspace(L,fac,n)
%
%   Generates a power series consisting of n points between which distance
%   grows or declines with a factor fac and starting at 0. The total length
%   of the series, i.e. sum(y) equals L.
%
%   with d = diff(y), then
%   The d(n)=d1*f^(n-1), and their sum is sum(d*f^(n-1)).
%
%
%   See also logspace linespace sinespace
%
%   TO 120118

if  nargin<3,
    n=100;
else
    n=round(double(n));
end

g=fac.^(0:n-1);

y=[0 cumsum(L*g/sum(g))]; y(end)=L;
