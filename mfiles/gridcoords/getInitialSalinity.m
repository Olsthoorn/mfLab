function C = getInitialSalinity(gr,Cmin,Cmax,xc,zc,sigmaX,sigmaZ)
%GETINITIALSALINITY generate sigmoid vertical salinity distribution
%
% Example:
%    C = getInitialSalinity(gr,Cmin,Cmax,xc,zc,sigmax,sigmaz)
%
% The salinity distribution is a sigmoid in z and x direction with
% Sigmax and Simagz respectively around point xc,yc
%
% See also: setEnvironmentalHead
%
% TO 120615 130416

% we assume 10 m fresh water diffusion depth blow the HLmeer lake
% hence sigmaZ.
if 0
    % Depth to compute erfc for conctration
    Depth = bsxfun(@minus,zc,XS(gr.ZM));

    % We will assume some horizontal transition of the sigmaZ around the
    % westeren shore of the Haarlem Lake.
    % Distance from HLlake west side
    sigmaZ = erfc((xc - gr.xm)/SigmaX) * SigmaZ;

    Cmean = 0.5*(Cmax-Cmin);

    %C=Cmin + Cmean* XS(erfc(- bsxfun(@rdivide, Depth,sigmaZ)));

    C=Cmin + Cmax* XS(erfc(- bsxfun(@rdivide, Depth,sigmaZ))-1);

else

    wLand = bsxfun(@times,ones(gr.Nz,1),0.5 + 0.5 * erf((gr.xm-xc)/sigmaX));
    wSea  = bsxfun(@times,ones(gr.Nz,1),0.5 + 0.5 * erf((xc-gr.xm)/sigmaX));

    cSea = Cmax * ones(gr.Nz,gr.Nx);
    cLand= bsxfun(@times,ones(1,gr.Nx), (Cmax - Cmin) *erf((zc - gr.zm(:))/sigmaZ));

    C = XS(wSea .* cSea + wLand .* cLand);

end


