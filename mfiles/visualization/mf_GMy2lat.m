function Lat=mf_GMy2lat(y)
%MF_GMY2LAT computes Latitude from googles y
%
% Example:
%   Lat=mf_GMy2lat(y);
% 
%   LAT=Northing
%
% SEE ALSO: mf_GMlat2y mf_GMlon2x mf_GMx2lon mf_GMLL2pix mf_GM2PIC mf_GM2PNG mf_GMpix2m
%
% TO 110506

phi=2*atan(exp(pi*(1-2*y)))-pi/2;

Lat=180/pi*phi;

