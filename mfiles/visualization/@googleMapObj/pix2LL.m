function [Lat Lon]=pix2LL(o,px,py)
% [Lat Lon] =GM.pix2LL: computes Lat Lon from local GM pix and tile coordinates
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

if nargin<5, [o.pxLL,o.pyLL]=ginput; end

w=2^-o.zoom;

y=w*(o.iy+(py+0.5)/255);
x=w*(o.ix+(px+0.5)/255);

Lat=GM.lat(y);
Lon=GM.lon(x);
