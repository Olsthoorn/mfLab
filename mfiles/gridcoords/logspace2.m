function [X,DX]=logspace2(x1,x2,dxMin,dxMult,dxMinBack)
%LOGSPACE2 generates a nice spacing based on end points
%
% Example:
%   [x,dx]=logspace2(x1,x2,dxMin,dsMult,dxMinBack)
%
%   Matlab's standard functions linspace and logspace facilitate generating
%   grids with a variable distances between grid lines.
%   logspace2 uses start and end locations and startdx and a dx multiplier.
%   dxMinBack is the starting value at x2 in case dx near x2 > than the
%   desired value of dxMinBack. Default is dxMinBack=dxMin
%
% See also: linespace logspace logspace2 gridObj gridTutorial
%
% TO 100214 100918 130701

%% Self test if nargin 0
if nargin==0
    selftest;
    return;
end

%% Assert proper input
if dxMult<=1, error('dxmult = %g, it must be > 1',dxmult); end

if nargin<5
    dxMinBack = abs(dxMin);
end

dxMin  = abs(dxMin);
dxMult = abs(dxMult);
L       =abs(x2-x1);

n1 = (log(L/2) -log(dxMin    ))/log(dxMult);
n2 = (log(L/2) -log(dxMinBack))/log(dxMult);

dx1=ones(1,max(1,fix(n1)+1))*dxMult;
dx2=ones(1,max(1,fix(n2)+1))*dxMult;

cdx1= cumsum(dxMin     * [1 cumprod(dx1)]);
cdx2= cumsum(dxMinBack * [1 cumprod(dx2)]);

X1 = x1+sign(x2-x1)*[0 cdx1];
X2 = x2+sign(x1-x2)*[0 cdx2];
if x2>x1
    X1 = X1(X1<=x1+L/2);
    X2 = X2(X2>=x2-L/2);
    if X2(end)-X1(end) > max(diff(X1));
        X = [X1 mean([x1 x2]) X2(end:-1:1)];
    else
        X = [X1 X2(end:-1:1)];
    end
else % x2<x1
    X2 = X2(X2<=x2+L/2);
    X1 = X1(X1>=x1-L/2);
    if X1(end)-X2(end) < max(diff(X2));
        X = [X1 mean([x1 x2]) X2(end:-1:1)];
    else
        X = [X1 X2(end:-1:1)];
    end
end
DX = abs(diff(X));

end

function selftest()
    x2 = 350;
    x1 = 35;
    dxMin = 0.5; dxMinBack = 0.1;
    dxMult= 1.25;
    display(logspace2(x1,x2,dxMin,dxMult));
    display(logspace2(x2,x1,dxMin,dxMult));
    display(logspace2(x1,x2,dxMin,dxMult,dxMinBack));
    display(logspace2(x2,x1,dxMin,dxMult,dxMinBack));
end
