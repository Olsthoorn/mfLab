function Lon=mf_GMx2lon(x)
%MF_GMX2LON computes googles x from Longitude
%
% Example:
%    Lon=mf_GMx2lon(x);
% 
%    lonN=Easting
%
% SEE ALSO: mf_GMx2lon mf_GMlat2y mf_GMy2lat mf_GMLL2pix mf_GM2PIC mf_GM2PNG mf_GMpix2m
%
% TO 110506

lam=pi*(2*x-1);

Lon=180*lam/pi;
