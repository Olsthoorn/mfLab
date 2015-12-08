function [X,DX]=logspace3(x1,x2,dxMin,dxMult,dxMax)
%LOGSPACE3 generates a logspace between x1 and x2 starting with dxMin using
% multiplyer dxMult and having dxMax as max dx size.
%
% Example:
%   [x,dx]=logspace2(x1,x2,dxMin,dsMult,dxMax)
%
%   Matlab's standard functions linspace and logspace facilitate generating
%   grids with a variable distances between grid lines.
%   logspace2 uses start and end locations and startdx and a dx multiplier.
%   dxMinBack is the starting value at x2 in case dx near x2 > than the
%   desired value of dxMinBack. Default is dxMinBack=dxMin
%
% See also: linespace logspace logspace2 logspace3 gridObj gridTutorial
%
% TO 100214 100918 131130

%% Self test if nargin 0
if nargin==0
    return;
end

%% Assert proper input
if dxMult<=1, error('dxmult = %g, it must be > 1',dxmult); end

if nargin<5
    dxMax = Inf;
end

dxMin  = abs(dxMin);  % dxMin must be > 0
dxMult = abs(dxMult); % dxMult must be >0
x1     = abs(x1);
x2     = abs(x2);
m      = min(x1,x2);
M      = max(x1,x2);
L       =M-m;  % x1 and x2 must be >0

n = fix(log(L)/log(dxMult))+1; % n is large enough

X = x1 + dxMin * [0 cumsum(dxMult.^(0:n))];
dx = diff(X);
dx= dx(dx<dxMax);

X  = x1+[0 cumsum(dx)]; X = [X(1:end-1)  X(end):dxMax:x2];

if X(end)<x2, X=[X x2]; end

DX = diff(X);

end
