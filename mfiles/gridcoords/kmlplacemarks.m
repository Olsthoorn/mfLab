function plmrk=kmlplacemarks(kmlfile,type,plmrk)
%KMLPLACEMARKS gives wgs-coodinates of all placemarks in GoogleEarth path embedded in kml file
%   useful after saving an entire folder from GE as kml files, which
%   contains all placemarks
%
%   pms=kmlplacemarks(kmlfname)
%   kmlfname = your kml file name as saved from GE
%   pms is a struct containing name, type and WGS coordiantes of placemarks
%   should work for any type of placemark in the kml file
%
%
% Example: (recepie)
%   1) launch GE, select add path, click a path
%   2) when finished save it as a kml file in a convenient directory
%      e.g. as "mypath.kml"
%   3) type
%   pms=kmlplacmarks('mypath');
%
% SEE ALSO  kmlpath2rd wgs3rd rd2wgs kmlpath
%
% TO 091215 110427 110501 111218 (LAT LON added)

types={'LineString','Polygon'};

if ~exist('type','var') || isempty(type)
    type=types{1};
    fprintf('Only one argument used, default type = <%s>!\n',type);
end

if isempty(strfind(kmlfile,'.kml')),
    kmlfile=[kmlfile '.kml'];
end

if iscell(type)
    for i=1:length(type)
        if ~exist('plmrk1','var')
            plmrk=kmlplacemarks(kmlfile,type{i});
        else
            plmrk=[plmrk;kmlplacemarks(kmlfile,type{i},plmrk)];
        end
    end
    return;
end

if ~strmatchi(type,types,'exact')
    for i=1:length(types), fprintf('%s\n',types{i}); end
    error('illegal type <%s>, use one of the just printed types!',type);
end    

%% read the entire kml file in as char for processing
fid=fopen(kmlfile,'r');
    if fid<1, error('Can''t open file %s\n',kmlfile); end
    s=fscanf(fid,'%c',[1,Inf]);
fclose(fid);

%% Seek location of tokens in this file/string
pm.Ipm     = strfind(s,'<Placemark>');
pm.IpmEnd  = strfind(s,'</Placemark>');
pm.Iname   = strfind(s,'<name>');
pm.InameEnd= strfind(s,'</name>');
pm.LType   = strfind(s,['<' type '>']);
pm.LTypeEnd= strfind(s,['</' type '>']);
pm.Coo     = strfind(s,'<coordinates>');
pm.CooEnd  = strfind(s,'</coordinates>');

%% Check heuristically to see that file is complete
if length(pm.Ipm) ~= length(pm.IpmEnd)
    error('kml file %s is incomplete, numer of <placemark> (%d) and of </placemark> (%d) are unequal!\n',...
        length(pm.Ipm),length(pm.IpmEnd));
end

%% Allocate memory for struct array to hold placemarks

plmrk2(length(pm.Ipm),1).type='dummy';

%% Fill this struct array
for i=1:length(plmrk2)
    
    % placemark name (a name is always present)
    j=find(pm.Iname>pm.Ipm(i) & pm.Iname<pm.IpmEnd(i)); % does placemark have a name token embedded?
    J=[pm.Iname(j) pm.InameEnd(j)];                       % must always be the case, error wil be thrown if not
    
    plmrk2(i).name=sscanf(s(J(1)+length('<name>'):J(2)-1),'%c',[1,Inf]); % scan the actual name of placemark

    % We will look for embedded type
    if any(pm.LType>pm.Ipm(i) & pm.LType<pm.IpmEnd(i))      % maybe it does not exist
        plmrk2(i).type=type;
    else
        plmrk2(i).type='none';
        fprintf('name: %s ---- type: %s\n',plmrk2(i).name,plmrk2(i).type);
    end
    
    % Get the coordinates (always present for each placemark)
    j=find(pm.Coo>pm.Ipm(i) & pm.Coo<pm.IpmEnd(i));
    J=[pm.Coo(j) pm.CooEnd(j)];
    
    ss=s(J(1)+length('<coordinates>'):J(2)-1);
    ss(ss==',')=' ';               % remove commas
    crds=sscanf(ss,'%f',[3,Inf])'; % scanf coords from string
    plmrk2(i).E=crds(:,1);         % get Easting
    plmrk2(i).N=crds(:,2);         % get Northing
    plmrk2(i).Z=crds(:,3);         % get elevation
end

%% remove all placemark of different types
I=strmatchi(type,{plmrk2.type},'exact');
if I~=0, plmrk2=plmrk2(I); end
if exist('plmrk','var'),
    plmrk=[plmrk;plmrk2];
else
    plmrk=plmrk2;
end

