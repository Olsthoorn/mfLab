function STCONC=salinityDistribution(varargin)
%SALINITYDISTRIBUTION -- generates a parameterized salinity distribution based on erf(x)
%
% generate parameterized start concentratrations and interface along the x-axes
% with sigmoid brackish water distribution using the error function
% sigma is the standard deviation defining the erfc function
%
% USAGE:
%       STCONC=salinityDistribution(gr,zeta,sigma,Cmax)
%
% zeta     = elevation of center of interface given by a set of points
% sigma    = standard deviation of interface thickness given by a set of points 
%
% TO 120510

[gr,   varargin] = getType(varargin,'gridObj',[]);
[zeta, varargin] = getNext(varargin,'struct',[]);
[sigma,varargin] = getNext(varargin,'double',[]);
[CMax ,~       ] = getNext(varargin,'double',[]);

if isempty(gr) || nargin<2
    STCONC = selftest();
    return;
end

method = 'spline';

Zeta = bsxfun(@times,interp1(zeta.x, zeta.y,gr.xm,method),ones(1,1,gr.Nz));
Sigma= bsxfun(@times,interp1(zeta.x,sigma,  gr.xm,method),ones(1,1,gr.Nz));
        
STCONC=CMax*erfc((gr.ZMlay-Zeta)./Sigma)/2;

figure; hold on; contourf(gr.xc,gr.zc,XS(STCONC),0:0.025:1,'edgeColor','none');
colorbar;

end

function STCONC = selftest()

CMax = 1;

xGr = -6500:100:20500;
yGr = [-0.5 0.5];
zGr = 0:-5:-250;
gr = gridObj(xGr,yGr,zGr);

zeta.x = [ -6  -4 -3 -2.5  -2     0    2    4   6     8  10  12   14  16  18   20] * 1000;
zeta.y = [100 100  0 -75  -100  -130 -125 -115 -105 -100 -95  -90 -90 -90 -90  -90];
sigma  = [  1   2  6  12    15    17   18   20  22.5  25  27.5  30 32.5 35  36   36.5];  

STCONC = salinityDistribution(gr,zeta,sigma,CMax);

end
