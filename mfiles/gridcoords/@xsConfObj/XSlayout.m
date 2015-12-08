function [values,Conf,mat]=XSlayout(o,basename,configSheetnm,materialSheetnm)
%MF_CONF: Gets grid values of the parameter parnam from specified zones
%   the zones are specified in the
%   worksheets "config" and "material" of Excel file [basename '.xls']
%   Used for convenient and efficient generation of cross sections.
%   See the examples that use it, e.g.
%        mflab>examples>swt_v4>freshkeeper>
%
% USAGE:
%   [values,Conf,mat]=mf_zone(basename,xGr,yGr,zGr,config,material,parnam)
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
% Material        Code	kh	    kv	red	green	blue	peff	rhofixed	rhodry	rhowet
% Plestocene_sand	P	30	    10	    1	1	0.5	0.35	2640	2066	2066
% Holocene_sand    	Z	1	    1	    1	0.75	0	0.35	2640	2066	2066
% Sludge	        S	0.02	0.02	0.2	0.8	0	0.5	1400	1200	1200
% Clay	            K	0.01	0.01	0	1	1	0.6	2100	1440	1440
% Peat	            V	0.04	0.04	0.5	0.5	0	0.6	1400	1160	1160
% Water	            W	1000000	1000000	0	0	0.8	1	1000	1000	1000
% Bottom sludge	    B	0.01	0.01	0.75	0.5	0	0.6	2100	1440	1440
% Air	            A	0	0	1	    1	1	1	0	0	
%
% to prevent missing some layers, make sure that the top and bottom
% elevations of the zones match with vertical grid lines. You can do that
% by merging the two sets of elevations when defining zGr. Just put the two
% elevation sets in one vector and use that in modelsize3(xGr,yGr,zGr),
% Modelsize3 guarantees correct merging.
%
% SEE ALSO: mf_zone mf_setwells mf_plotConf
%
% TO 110319
basename = 'Vechtbaggeren';
configSheetnm = 'Configuration';
materialSheetnm = 'Material';

%% Get material
[mat.header,mat.values,txtHdr,txtValues]=getExcelData(basename,materialSheetnm,'Hor');

mat.name =txtValues(:,strmatchi('Material',txtHdr));
mat.code =cell2mat(txtValues(:,strmatchi('Code'    ,txtHdr)));

%% Get zones

[zone.name,zone.values]=getExcelData(basename,configSheetnm,'Hor');

zone.Left  = zone.values(1,2); if isnan(zone.Left), zone.Left = 0; end
zone.width = zone.values(1,3:end);
zone.head  = zone.values(2,3:end);
zone.top   = zone.values(3,3:end);
zone.xL    = zone.Left +[0 cumsum(zone.width(1:end-1))];
zone.xR    = zone.Left +   cumsum(zone.width);
zone.bot   = zone.values(4:end,2);

[~,~,~,zone.code]=getExcelData(basename,configSheetnm,'Ver');
zone.code = cell2mat(zone.code(5:end,2:end));

%% replace zone codes by direct index into material list
Idx = (1:max(mat.code))'*[1 NaN];
Idx(mat.code(:,1),2)=(1:length(mat.code))';

zone.imat = reshape(Idx(zone.code,2),size(zone.code));

%% reproduce the grid

% Get zone index in full grid
Ix=ones(size(xm)); for ix=1:length(Conf.xL),     Ix(xm>Conf.xL( ix))=ix; end
Iz=ones(size(zm)); for iz=length(Conf.bot):-1:1, Iz(zm>Conf.bot(iz))=iz; end

%% get the parameter column
ipar = strmatchi(parnam,mat.header); % which parameter column:

%% fill full grid with parameter value

% any value outside the Conf layout will be a NaN
for ix=1:length(Conf.names)
    I=find(xm>Conf.xL(ix) & xm<Conf.xR(ix)); % ix indices model in this conf zone
    zConf = [Conf.Top(ix); Conf.bot(:,ix)];
    for i=1:length(I)
        values(:,I(i))=Transfergrids(zConf,ConfValues,gr.Z(1,I(i),:),code);
    end
   %     values(:,I,J)=mat.values(mat.idx(iz,ix),ipar);
   % end
end


% any value outside the Conf layout will be a NaN
for ix=1:length(Conf.names)
    I=find(xm>Conf.xL(ix) & xm<Conf.xR(ix));
    Conf.Head(I)=Conf.head(ix);
    Conf.Top(I) =Conf.top( ix);
    for iz=length(Conf.bot):-1:1
        if iz==1, J=find(zm>Conf.bot(1)  & zm<Conf.top(ix));
        else      J=find(zm>Conf.bot(iz) & zm<Conf.bot(iz-1));
        end
        values(:,I,J)=mat.values(mat.idx(iz,ix),ipar);
    end
end

for ix=1:gr.Nx
    values(1,ix,:)= gridsTransfer([Conf.top(iC); Conf.bot(ix)],Conf.param(),gr.zgr(1,ix,:),'k');
end

    
