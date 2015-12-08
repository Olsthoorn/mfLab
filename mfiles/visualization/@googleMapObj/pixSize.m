function [wLat wLon]=pixSize(o)
% [wLat wLon] = pixSize(o); Computes size of pixel in m in tile of center
%
% TO 110503 121008

R2 = 6371007.2;  % The Earth's authentic radius (Geodetic Union);

w=2^-o.zoom;

wLat =   pi*R2*w/256;
wLon = 2*pi*R2*w/256*cos(pi*o.center(1)/180);
