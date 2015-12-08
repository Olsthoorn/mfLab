function [dx,dy,r] = dist(Lat1,Lon1,Lat2,Lon2)
%DIST computes distance in m between two googleMap points
%
% Example
%    [dx,dy,r] = dist(lat1,lon1,lat2,lon2);  computes distance in m between two googleMap points
%    [r      ] = dist(lat1,lon1,lat2,lon2) compute radial distance if nargout = 1
%
% TO 110501 121008

R2 = 6371007.2;  % The Earth's authentic radius (Geodetic Union);

if nargin<2
    [Lon2 Lat2]=ginput;
end
gp1 = googlePointObj(Lat1,Lon1(1:max(1,numel(Lat1))));
gp2 = googlePointObj(Lat2,Lon2(1,max(1,numel(Lat2))));

if nargout>1
    [dx,dy,r] = gp1.dist(gp2);
else
    r         = gp1.dist(gp2);
end
