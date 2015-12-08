function conc = transitionZone(o,cmin,cmax,z50,sigma)
%% conc = transitionZone(o,cTop,cBot,z50,sigma)
%  generates a smooth salinity array based on the error function
%  cmin is minimum conc
%  cmax is maximum conc
%  z50  is elevation of 50% plane
%  sigma is standard deviation used in the error function, a good masure of
%  the half thickness of the transition zone.
%  z50 and sigma may be sclars or 2D arrays of size [gr.Ny, gr.Nx]
%
%  TO DO: verify TO 120807
%
% TO 120807 130408

z50   = bsxfun(@times,ones(1,1,o.Nlay)*z50  ,ones(o.Ny,o.Nx));
sigma = bsxfun(@times,sigma,ones(o.Ny,o.Nx));
sigma = bsxfun(@times,sigma,ones(o.size));

conc = 0.5*(cmin + cmax) + 0.5*(cmax - cmin) * erf((z50-o.ZM)./sigma);
