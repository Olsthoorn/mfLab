function [X,Y]=mf_kmlpath2rd(fname)

%   KMLPATH: Give wgs-coodinates  of GE path in kml file
%
%   [E,N]=kmlpath(kmlfname)
%   kmlfname = fname given by you to your savd as kmfl file Google Earth path
%
%   LL = kmlpath(kmlfname)
%   with one output argument the coordinates are retured as matrix [Lat Lon]=[N E]
%
%   wgs coordinates are Easting and Northing in that order
%
%EXAMPLE
%   launch GE, select add path, click a path
%   when finished save it as a kml file in a convenient directory
%   e.g. as "mypath.kml"
%   type
%   [N,N]=mf_kmlpath('mypath');
%
%SEE ALSO  mf_kmlpath2rd wgs3rd rd2wgs
%
% TO 091215 110427 110501 (LAT LON added)

[P N E]=fileparts(fname);
 if isempty(E), fname=[fname '.kml']; end
 
fp=fopen(fname,'r');
if fp<0, error('Can''t open file <<%s>>',fname); end

s=fgetl(fp);
while s~=-1
    if findstr('<coordinates>',s)
        s=fgetl(fp);                  % coordinates are on next line
        i=findstr('</co',s); 
        s(i:end)='';                  % remove trail
        s(s==',')=' ';                % remove commas
        crds=sscanf(s,'%f',[3,Inf])'; % scanf coords from string
        E=crds(:,1);                  % get Easting
        N=crds(:,2);                  % get Northing
        break;
    end
    s=fgetl(fp);
end

[X,Y]=wgs2rd(E,N);

if nargout==1, X=[X Y]; end 

fclose(fp);
