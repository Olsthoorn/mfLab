function [Lat Lon]=pix2LL(ix,iy,px,py,zoom)
%PIX2LL convert PIX coordinates to Lat Lon
% [Lat Lon] =pix2LL(ix,iy,px,py,zoom): computes Lat Lon from local google point coordinates
%
% USAGE:
%   [Lat Lon]=GM.pix2LL(px,py); 
%
%   ix=GM x-tile number (0 .. (2^n)-1
%   iy=GM y-yile number (0 .. (2^n)-1
%   px=GM within tile coordinate (0-255)
%   py=GM within tile coordinate (0-22)
%
% SEE ALSO: GM.LL2pix
%
% TO 110501 121008

w=2^-zoom;

y=w*(iy+(py+0.5)/255);
x=w*(ix+(px+0.5)/255);

phi=2*atan(exp(pi*(1-2*y)))-pi/2;
Lat=180/pi*phi;

lam=pi*(2*x-1);
Lon=180*lam/pi;

