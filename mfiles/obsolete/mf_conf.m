function [values,Conf,mat]=mf_conf(basename,xGr,yGr,zGr,config,material,parnam)
% MF_CONF: Gets grid values of the parameter parnam from specified zones
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

%% Get material
[Num,Txt,Material]=xlsread(basename,material,'','basic');

mat.header=Material(1,3:end);
mat.names =Material(2:end,1);
mat.codes =Material(2:end,2);
mat.values=NaN(size(Material(2:end,3:end)));
for j=1:size(mat.values,1)
    mat.values(j,:)=horzcat(Material{j+1,3:end});
end

%% Get zones
[Num,Txt,Raw]=xlsread(basename,config,'','basic');
Conf.names= Raw(1,3:end);
Conf.xL = horzcat(Raw{2,2:end-1}); Conf.xR=horzcat(Raw{2,3:end});
Conf.head = horzcat(Raw{3,3:end});
Conf.top  = horzcat(Raw{4,3:end});
Conf.bot  = vertcat(Raw{5:end,2});
Conf.mzone= Raw(5:end,3:end);

mcodes =unique(vertcat(Conf.mzone(:)));

mat.idx=NaN(size(Conf.mzone));

for i=1:length(mcodes),
    mat.idx(strmatchi(mcodes{i},Conf.mzone))=strmatchi(mcodes{i},mat.codes);
end

%% reproduce the grid
[xGr,yGr,zGr,xm,ym,zm,Dx,Dy,Dz,Nx,Ny,Nz]=modelsize3(xGr,yGr,zGr);

% Get zone index in full grid
Ix=ones(size(xm)); for ix=1:length(Conf.xL),     Ix(xm>Conf.xL( ix))=ix; end
Iz=ones(size(zm)); for iz=length(Conf.bot):-1:1, Iz(zm>Conf.bot(iz))=iz; end

%% get the parameter column
ipar = strmatchi(parnam,mat.header); % which parameter column:

%% fill full grid with parameter value
values=zeros(Ny,Nx,Nz);
conf.Top =NaN(Nx);
conf.Head=NaN(Nx);

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
