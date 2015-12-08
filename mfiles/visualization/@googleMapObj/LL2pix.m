function gp = LL2pix(o,Lat,Lon)
% GMLL2PIX: translates lat lon to pix coordinates on local GM tile
%
% USAGE:
%   gp=GM.LL2pix(Lat,Lon)
%
%   gp=googlePointObj
%
%   input:N= LAT E=LON and the used o.zoom level 0-21
%   px,py are pixel coordinates from the LAT LON pont within the GM tile
%
%   ix,iy is the tile coordinates 0 ..(2^n)-1 level/position in Goolgle Maps tile system
%
% TO 110501 121008

gp.w=2^-o.zoom;

y=o.y(Lat);
x=o.x(Lon);

gp.ix=fix(x/gp.w);
gp.iy=fix(y/gp.w);

gp.px =255*(x/w-gp.ix);
gp.py =255*(y/w-gp.iy);