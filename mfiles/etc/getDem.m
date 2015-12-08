function Z=getDem(tifFName)
%GETDEM gets dem from USGS tiffile specific for MS reservoir
%
% USAGE:
%    getDem(tifFName) --- 
%
% project model grid
%
% Copyright Hessel Winsemius 070501

[dem,R,BBOX]=geotiffread(tifFName);

latlim = BBOX(:,2);
lonlim = BBOX(:,1);

mstruct = defaultm('utm');
mstruct.zone = '23L';
wgs84 = almanac('earth','wgs84');
format long g
wgs84(1) = wgs84(1)*1000;
mstruct.geoid = wgs84;
mstruct = defaultm(mstruct);

[xlim ylim] = mfwdtran(mstruct,latlim,lonlim);

xax = 239000:  100:246000;
yax = 8243000:-100:8235900;
[xi,yi] = meshgrid(xax,yax);
[lati,loni] = minvtran(mstruct,xi,yi);

lat2 = linspace(-15.8779-0.00083333333333/2,-15.9421+0.00083333333333/2,77);
lon2 = linspace(-47.4371+0.00083333333333/2,-47.3721-0.00083333333333/2,78);

[loni2,lati2] = meshgrid(lon2,lat2);

Z = interp2(loni2,lati2,dem,loni,lati);

%figure
%imagesc(xax,yax,Z)
%figure
%surf(double(xax), double(yax), double(flipud(Z)))
%set(gca,'YDir','normal')
