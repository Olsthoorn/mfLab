function Lon_=lon(~,x)

% lon_ = GM.long(x): Computes Google's x from Longitude
%
% USAGE:
%   Lon=GM.x;
% 
%   Lon=Easting
%
% SEE ALSO: GM.lat GM.y GM.x mf_GMLL2pix mf_GM2PIC mf_GM2PNG mf_GMpix2m
%
% TO 110506 121007

lam=pi*(2*x-1);

Lon_=180*lam/pi;
