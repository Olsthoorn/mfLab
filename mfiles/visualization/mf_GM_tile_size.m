function [dx,dy]=mf_GM_tile_size(zoom,Lat,varargin)
%MF_GM_TILE_SIZE computes size of pixel in m according to zoom level
%
% Example:
%   [dx,dy]=mf_GM_tile_size(zoom,Northing);
%   [dE,dN]=mf_GM_tile_size(zoom,Northing,'EN');
%
%   [dx,dy] = tile size in m
%   [dE,dN] = tile size in E and N (degrees)
%
% SEE ALSO: mf_GMLL2pix mf_GM2PIC mf_GM2PNG mf_GMpix2m
%
% TO 110503

R2 = 6371007.2;  % The Earth's authentic radius (Geodetic Union);
n=2^zoom;
w=1/n;

iy=fix(mf_GMlat2y(Lat)/w);
y1=w*iy;
y2=w*(iy+1);

N1=mf_GMy2lat(y1);
N2=mf_GMy2lat(y2);

dy=R2*pi/180*(N1-N2);

dx=w*2*pi*R2*cos(pi*Lat/180);

if nargin>2
    dy=N1-N2;
    dx=w*360;
end