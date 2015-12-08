function [Values,Conf,material]=materialObj(basename,xGr,yGr,zGr,config,mat,parnam)
%MATERIALOBJD replaced by xyConfObj
%
% Under constructin unclear whether this is a useful object TO 120613
%
% MF_ZONE: Gets grid values of parameter parname from zones specified in the worksheet
%    with name config and material of Excel file [basename '.xls']
%
% USAGE:
%   [values,Conf,material]=mf_zone(basename,xGr,yGr,zGr,config,material,parnam)
%
% zGr must be one 3D vector zGR(1,1,:).
%
% Layout of config worksheet:
%  line 1 afer Base(m) are arbitrary zone names (arbitrary)
%  line 2 first col is label, then follow left x of left-most zone, and
%  right x of all zones
%  line 3: label, blank, top of all zones
%  line 4: label, blank, head_of all zones
%  line 5: layer nr, layer bottom, character code id for all zones
%  line 6: same
%  line 7: same etc
%
%Example:
%Layer	Base(m)	KWA	UITH	AMST	EBDK	EBDK	GMIJ	VINKV	ZGAT	VNKV	EAST
%xZone	-8000	-5200	-4500	-4000	-2700	-1800	1400	3000	3500	4500	7000
%head_top	-5.0	-3.0	-1	-6.0	-6.0	-6.5	-2	-2	-2	-2.0
%Top_(m)	-4.5	-2.5	-1	-5.5	-5.5	-6.0	-2	-2	-2	-1.5
% 1	  -8	V	V	V	V	V	V	V	V	V	V
% 2	 -40	P	P	P	P	P	P	P	P	P	P
% 3	 -50	K	K	K	K	P	P	P	P	P	P
% 4	-250	P	P	P	P	P	P	P	P	P	P
% 5	-260	K	K	K	K	K	K	K	K	K	K
% 6	-350	P	P	P	P	P	P	P	P	P	P
%
% The layout of the material worksheet
% line 1; label, label names of all material properties (arbitrary)
% line 2: full name, character code, values
% line 3: same
% line 4: ame etc
%
% Example:
% Material        Code	kh	kv	red	green	blue	peff	rhofixed	rhodry	rhowet
% Plestocene_sand	P	30	10	1	1	0.5	0.35	2640	2066	2066
% Holocene_sand    	Z	1	1	1	0.75	0	0.35	2640	2066	2066
% Sludge	            S	0.02	0.02	0.2	0.8	0	0.5	1400	1200	1200
% Clay	            K	0.01	0.01	0	1	1	0.6	2100	1440	1440
% Peat	            V	0.04	0.04	0.5	0.5	0	0.6	1400	1160	1160
% Water	            W	1000000	1000000	0	0	0.8	1	1000	1000	1000
% Bottom sludge	    B	0.01	0.01	0.75	0.5	0	0.6	2100	1440	1440
% Air	            A	0	0	1	1	1	1	0	0	
%
% To prevent missing some layers, make sure that the top and bottom
% elevations of the zones match with vertical grid lines. You can do that
% by merging the two sets of elevations when defining zGr. Just put the two
% elevation sets in one vector and use that in modelsize3(xGr,yGr,zGr),
% Modelsize3 guarantees correct merging.
%
% SEE ALSO: mf_conf, mf_setwells mf_plotConf
%
% TO 110319

%% Get material
warning('off','all');
[~,~,Material]=xlsread(basename,mat,'','basic');
warning('on','mf_zone:XLSREAD:basic');

material.header=Material(1,3:end);
material.names =Material(2:end,1);
material.codes =Material(2:end,2);
material.values=NaN(size(Material(2:end,3:end)));
for j=1:size(material.values,1)
    material.values(j,:)=horzcat(Material{j+1,3:end});
end

%% Get zones
warning('off','all');
[~,~,Raw]=xlsread(basename,config,'','basic');
warning('on','mf_zone:XLSREAD:basic');

dir=lower(Raw{2,1}(1));
dL=[dir 'L'];
dR=[dir 'R'];

Conf.names= Raw(1,3:end);
if dir=='X'
    Conf.(dL) = horzcat(Raw{2,2:end-1}); Conf.(dR)=horzcat(Raw{2,3:end});
else
    Conf.(dL) = horzcat(Raw{2,2:end-1}); Conf.(dR)=horzcat(Raw{2,3:end});
end

if dir~='x'; error('Dir must be ''x'', other direction not (yet) implemented!'); end

Conf.head = horzcat(Raw{3,3:end});
Conf.top  = horzcat(Raw{4,3:end});
Conf.bot  = vertcat(Raw{5:end,2});
Conf.mzone= Raw(5:end,3:end);

mcodes =unique(vertcat(Conf.mzone(:)));

material.idx=NaN(size(Conf.mzone));

for i=1:length(mcodes),
    material.idx(strmatchi(mcodes{i},Conf.mzone))=strmatchi(mcodes{i},material.codes);
end

%% Join with model grid (merging)
%% if zGr for all columns is equal, then we can merge the elevations
%  of the zones with those of the zGr
%  we will assume this here
Merged_Z  = [zGr(:); Conf.bot; Conf.top'];
Merged_gr = gridObj(xGr,yGr,Merged_Z);
Merged_ZM =XS(Merged_gr.ZMlay);

% Iz_Conf is the vertical   Conf     zone index for each model cell
% Ix_Conf is the horizontal Confzone index for each model cell

% first compute the Iz_Conf zone for all cells in the cross sections

% ZC are the layer elevations of the Conf zones / material zones
ZConf  = NaN(size(Conf.mzone,1)+1,size(Conf.mzone,2));
ZConf(2:1+size(Conf.mzone,1),:) = Conf.bot*ones(size(Conf.top));
ZConf(1,:) = Conf.top;
Iz   = (1:size(ZConf,1))';  % to interpolate to get the cell index

IzConf = NaN(size(Merged_ZM));
for i=1:length(Conf.top)
    Ix=find(Merged_gr.xm>Conf.xL(i) & Merged_gr.xm<Conf.xR(i));
    % This works even if each ZM vertical column is different,
    % and even if also the ZC of each conf zone is different !!
    IzConf(:,Ix)=floor(interp1(ZConf(:,i),Iz,Merged_ZM(:,Ix)));
end

% Get zone index in full grid
XConf=[Conf.xL  Conf.xR(end)];
Ix=1+(0:size(Conf.xR,2));
Merged_XM=XS(Merged_gr.XMlay);
IxConf=floor(interp1(XConf,Ix,Merged_XM));

%% get the material values pertaining te the specified parameter parnam
%  and put is in an array the size of the Conf zones
ipar = strmatchi(parnam,material.header); % which parameter column:
parvals    = material.values(:,ipar);
zoneValues = parvals(material.idx);

% Tackling model cells out of the Conf zone boundaries
% Dummy zone to bottom and end of zoneValues to allow a value
% for the mentioned model cells
zoneValues(end+1,:)=0;  % better than NaN at least for kh and kv
zoneValues(:,end+1)=0;  % better than NaN

% Make sure the out of reference NaNs in Ix_Conf and Iz_Conf point into
% the dummy zones to get their dummy value for this material
IzConf(isnan(IzConf))=size(zoneValues,1);
IxConf(isnan(IxConf))=size(zoneValues,2);

%% The array has the size of the merged ZM array !
Idx=size(zoneValues,1)*(IxConf(:)-1)+IzConf(:);
Merged_Values=zoneValues(Idx);
Merged_Values=reshape(Merged_Values,size(IzConf));
Merged_DZ = XS(Merged_gr.DZlay);
Merged_ZM = XS(Merged_gr.ZMlay);

%% new grid for original matrices
Orig_gr = gridObj(xGr,yGr,zGr);
Orig_Z  = XS(Orig_gr.Zlay);
Orig_DZ = XS(Orig_gr.DZlay);
Orig_Iz =(1:size(Orig_Z,1))';

Merged_IZ =floor(interp1(Orig_Z(:,1),Orig_Iz,Merged_ZM)); % assuming all Oriz_Z are the same (for sake of interpolation)

%% Now, finally, get back to the array of the mmodel

% Filling the original arrays additively line by line using these indices
Values=XS(NaN(Orig_gr.size()));

kv=strcmpi('kv',material.header{ipar});
for iz=1:size(Merged_IZ,1)
    if kv
        values=Merged_DZ(iz,:)./Merged_Values(iz,:);
    else
        values=Merged_Values(iz,:).*Merged_DZ(iz,:);
    end
    J = find(~isnan(Merged_IZ(iz,:)));
    if ~isempty(J)
        Idx=size(Values,1)*(J-1)+Merged_IZ(iz,J);
        Values(Idx) = values(J);
    end
%     for j=1:size(values,2)
%         Values(Merged_IZ(iz,j),j)=values(j);
%     end
end

Values=Values./Orig_DZ;

if kv
    Values=1./Values;
end

Values=XS(Values);

