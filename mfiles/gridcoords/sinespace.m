function [x,dx]=sinespace(x1,x2,N,alfa1,alfa2)
%SINESPACE generates a nice spacing based on end points
%
% Example:
%   [x,dx]=sinespace(x1,x2,N,alfa1 [,alfa2])
%
%   Matlab's standard functions linspace and logspace facilitate generating
%   grids with a variable distances between grid lines.
%   Sinespace uses the sine function to generate grid distances using alfa1 and
%   optionally alfa2, in radians, to condition the sine function on the end
%   points x1 and x2, while N is the number of grid lines including x1 and x2
%   sinespace(x1,x2,N,alfa) is equivalent to sinespace(x1,x2,N, 0,alfa)
%
% See also: linespace logspace logspace2 gridObj gridTutorial
%
% TO 100214 100918

if nargin==4, alfa2=alfa1; alfa1=0; end

L=abs(x2-x1);

u=abs(sin(alfa1+(alfa2-alfa1)*(1:N-0.5)/N)); % compute alfa for dx

dx=L*u/sum(u);     % match distances to total distance, dx always positive

if x2>x1
    x=[x1 x1+cumsum(dx)]; x(end)=x2;
else
    x=[x1 x1-cumsum(dx)]; x(end)=x2;
end
