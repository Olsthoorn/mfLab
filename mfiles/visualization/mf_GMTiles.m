function mf_GMTiles(lat,lon,zoomlev)
%MF_GMTILES translates lat lon to pix coordinates on local GM tile
%
% ToDO: check description above
%
% EXAMPLE:
%   [xp,yp,ix,iy]=GMLL2pix(lat,lon,zoomlev)
%
%   lat lon as usuze zoom leve 0-21
%   xp,yp pixel coordinates in local tile whose
%   ix,iy tile position in Goolgle Maps tile system
%   zoomlev = 0-21 according to google
%
% DETAILS:
%    Thanks to the math given by
%    Markus Loecher, Sense Networks, markus@sensenetworks.com Dec 3, 2010
%    see his pdf on the web.
%
% SEE ALSO: mf_GM2PIC mf_GMpix2LL
%
% TO 110501

zoom=10;
zoomfac=2^(zoom-1);

lon=179.99; LON=pi*lon/180;
lat= -0.01; LAT=pi*lat/180;

% X0,Y0 are google map tile coordinates with 0<X<1 and 0<Y<1 in tile 0 (all world)
Xaux=LON/180;
Yaux=1/(2*pi)*log((1+sin(LAT))/(1-sin(LAT)));

% X, Y are coordanates for te local tile such that the integer part is the
% tile number and the fractonal part the in-tile relative coordinate.
%
X=zoomfac*(Xaux+1);
Y=zoomfac*(1-Yaux);

% Tile number in google tile system
ix=floor(X);
iy=floor(Y);

% x,y are google map coordinats within local tile at zoom level n, also
% between 0 and 1
x=256*(X-ix);
y=256*(Y-iy);

end

function [LAT LON]=GMpix2LL(x,y,ix,iy,zoomlev)
% GMpix2LL: computes lat long from local pix coordinates given GM tile and
% zome level.
%
% USAGE:
%   [LAT LON]=GMpix2LL(x,y,ix,iy,zoomlev)
% 
%   ix=GM x-tile number
%   iy=GM y-yile number
%   may be obtained usgin [x,y,ix,iy]=GMLL2pix(Lat,Lon,zoomlev);
%
% SEE ALSO: mf_GMLL2pix mf_GM2PIC

X=ix+x/256;
Y=iy+y/256;

Xaux=X/zoomfac-1;
Yaux=1-Y/zoomfac;

lon=X/zoomfac-1;
lat= sin( (exp(2*pi*Yaux)-1) / (exp(2*pi*Yaux)+1) )+pi;

LON=180*lon;
LAT=180*lat;

end
 