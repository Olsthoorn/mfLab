function x_=x(~,Lon)

% x=GM.x(Lon): Computes googles x from Longitude
%
% USAGE:
%   x=mf_GMlon2x(Lon);
% 
%   lonN=Easting
%
% SEE ALSO: mf_GMLL2pix mf_GM2PIC mf_GM2PNG mf_GMpix2m
%
% TO 110506

lam=pi/180*Lon;

x_=(lam/pi+1)/2;

