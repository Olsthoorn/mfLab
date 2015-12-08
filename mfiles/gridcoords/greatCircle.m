function D=greatCircle(lat1,lon1,lat2,lon2,unit)
%GREATCIRCLE gets distance given Lat Long for two points
%
% Example:
%    D=greatCircle(lat1,lon1,lat2,lon2,unit);
%
% where unit is oneof 'k' (km) 'm' (miles) 'n' (nautical miles)
% notice: 'k' is also the default if unit is not specified.
%
% after: http://www.igeocode.com/developer/reference/great-circle-distance
%
% See also: rd2wgs wgs2rd googleMap
%
% TO 090101

if nargin<4
    error('%s: 4 inputs required, see help %s',mfilename,mfilename);
end

if nargin<5
    unit = 'k';
end

deg2rad = inline('d/180*pi','d');
radius = 6371.0090667; % earth mean radius by IUGG

dlon = lon1 - lon2;

distance = acos( sin(deg2rad(lat1)) * sin(deg2rad(lat2)) + ...
    cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * cos(deg2rad(dlon))) * radius; 

  switch unit
      case 'm', D=distance * 0.621371192;
      case 'n', D=distance * 0.539956803;          
      otherwise
          D=distance;
  end
end
