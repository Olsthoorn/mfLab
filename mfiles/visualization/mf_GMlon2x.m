function x=mf_GMlon2x(Lon)
%MF_GMLON2X computes googles x from Longitude
%
% Example:
%   x=mf_GMlon2x(Lon);
% 
%   lonN=Easting
%
% SEE ALSO: mf_GMLL2pix mf_GM2PIC mf_GM2PNG mf_GMpix2m
%
% TO 110506

lam=pi/180*Lon;

x=(lam/pi+1)/2;

