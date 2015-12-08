function [E_lon,N_lat,xUTM,yUTM]=kmlpath(fname)
%KMLPATH gives wgs-coodinates of GoogleEarth path embedded in kml file
%
% Example:
%     [E_lon,N_lat,xUTM,yUTM]=kmlpath(kmlfname)
%
%  where
%     kmlfname = fname given by you to your savd as kmfl file Google Earth path
%
%   LL = kmlpath(kmlfname)
%   with one output argument the coordinates are retured as matrix [Lat Lon]=[N E]
%
%   wgs coordinates are Easting and Northing in that order
%
% Example: (recepie)
%   1) launch GE, select add path, click a path
%   2) when finished save it as a kml file in a convenient directory
%      e.g. as "mypath.kml"
%   3) type
%     [N,E]=kmlpath('mypath');
%
% See also: kmlpath2rd wgs3rd rd2wgs
%
% TO 091215 110427 110501 (LAT LON added)

if isempty(strfind(fname,'.kml')), fname=[fname '.kml']; end
 
fp=fopen(fname,'r'); if fp<0, error('Can''t open file <<%s>>',fname); end


s=fscanf(fp,'%c',Inf); k=0;

n = regexp(s,'</?coordinates>');

crds=s(n(1)+length('<coordinates>'):n(2)-1);
crds(crds==',')=' ';
crds=sscanf(crds,'%f',[3,Inf])';
E_lon=crds(:,1);                  % get Easting
N_lat=crds(:,2);                  % get Northing

[xUTM,yUTM] = wgs2utm(N_lat,E_lon);

if nargout<2, E_lon=[E_lon,N_lat]; end  % [LAT LON] = [N E]

fclose(fp);
