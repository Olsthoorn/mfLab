function Lat_=lat(~,y)

% MF_GMY2LAT: Computes Latitude from googles y
%
% USAGE:
%   Lat=GM.lat(y);
% 
%   LAT=Northing
%
% SEE ALSO: googleMapObj.y googleMapObj.x googleMapObj.lon googleMapObj.pxLL googleMapObj
%
% TO 110506 1210-7

phi=2*atan(exp(pi*(1-2*y)))-pi/2;

Lat_=180/pi*phi;

