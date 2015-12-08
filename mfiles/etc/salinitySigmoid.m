function C = salinitySigmoid(gr,Cmin,Cmax,zc,sigmaZ)
%salinitySigmoid --- generate sigmoid vertical salinity distribution
%
% Example:
%    C = salinitySigmoid(gr,Cmin,Cmax,zc,sigmaz)
%
%  gr     = gridObj
%  Cmin   = min conc
%  Cmax   = max conc
%  xc     = coordinates or scalar
%  zc     = elevation of c_avg at xc
%  sigmaZ = std of sigmoid at xc
%
% TO 120615 130416

if isscalar(zc)
    zc     = ones(1,gr.Nx)*zc;
end
if isscalar(sigmaZ)
    sigmaZ = ones(1,gr.Nx)*sigmaZ;
end

% Depth to compute erfc for conctration
Depth = bsxfun(@minus,zc,XS(gr.ZM));

C= Cmin + (Cmax - Cmin)/2 * XS(erfc(- bsxfun(@rdivide, Depth,sigmaZ)));



