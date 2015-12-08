function LLstruct=kmlFolder(fname)
%KMLFOLDER Get path names and wgs coordinates of paths in kmlFolder
%
% Example:
%     LLstruct =kmlFolder(kmlfname)
%
%  where
%     kmlfname = fname given by you to your saved as kmfl folder file Google Earth path
%     LLstruct is a struct array with the following fields
%     LLstruct.name (path name)
%     LLstruct.E (wgs84 longitude (easting))
%     LLstruct.N (wgs84 latitude  (northing))
%     LLstruct.x (x UTM)
%     LLstruct.y (y UTM)
%
%   LL = kmlFolder(kmlfname)
%   with one output argument the coordinates are retunred as matrix [Lat Lon]=[N E]
%
%   wgs coordinates are Easting and Northing in that order
%
% Example: (recepie)
%   1) launch GE, select add paths in a folder in GE
%   2) add paths each with its own name
%   2) when finished save it as a kml file in a convenient directory
%      e.g. as "profiles.kml"
%   3) type
%     LLstruc=kmlFolder('mypath');
%
% See also: kmlpath kmlPath2UTM wgs2utm utm2wgs kmlpath2rd wgs3rd rd2wgs
%
% TO 091215 110427 110501 130926

if isempty(strfind(fname,'.kml')), fname=[fname '.kml']; end
 
fp=fopen(fname,'r'); if fp<0, error('Can''t open file <<%s>>',fname); end
fseek(fp,0,-1);

s=fscanf(fp,'%c',Inf);

if ~regexp(s,'<folder>','once')
    error('%s: File <<%s>> is not a kmlFolder file',mfilename,fname);
end

firstPlacemark = regexp(s,'<Placemark>','once');

nm1  = regexp(s,'<name>') + length('<name>'); % skip folder name
nm2  = regexp(s,'</name>')-1;

nm1 = nm1(nm1>firstPlacemark);
nm2 = nm2(nm2>firstPlacemark);

nCrd1  = regexp(s,'<coordinates>') + length('<coordinates>'); 
nCrd2  = regexp(s,'</coordinates>')-1;

nCrd1 = nCrd1(nCrd1>firstPlacemark);
nCrd2 = nCrd2(nCrd2>firstPlacemark);

for i=numel(nm1):-1:1
    LLstruct(i).name = s(nm1(i):nm2(i));
    crds= s(nCrd1(i):nCrd2(i));
    crds(crds==',')=' ';
    crds=sscanf(crds,'%f',[3,Inf])';
    LLstruct(i).E=crds(:,1);                  % get Easting
    LLstruct(i).N=crds(:,2);                  % get Northing
    [LLstruct(i).X,LLstruct(i).Y] = wgs2utm(LLstruct(i).N,LLstruct(i).E);
end

fclose(fp);
