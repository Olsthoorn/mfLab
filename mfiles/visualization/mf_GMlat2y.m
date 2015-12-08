function y=mf_GMlat2y(Lat)
%MF_GMLAT2Y computes googles Y from Laititute
%
% Examle:
%   y=mf_GMlat2y(Lat);
% 
%   LAT=Northing
%
% SEE ALSO: mf_GMLL2pix mf_GM2PIC mf_GM2PNG mf_GMpix2m
%
% TO 110506

phi=pi*Lat/180;

y=0.5*(1-log((1+sin(phi))./(1-sin(phi)))/(2*pi));

