%% Input for analytic and numeric modeling of varying average head in cross sections

% Dutch top system consisting of a shallow Holocene aquifer underlain by a
% semi-pervious layer on top of a larger regional aquifer. The shallow
% aquiifer is intersected by may dithces in which the surface water level
% is maintained to manage the groundwater system.
% We simulate two parallel ditches in this shallow top system underlain by
% the confining bed and the regional system. THe ditches may cut through
% the confining layer.
% Input is recharge and seepage from below through the confining bed from
% the regional aquifer. Extraction is evapotranspiration and leakage to the
% ditches. These dithches may at times also infiltrate.
% The specific yield will be made dependent on the the water table
% elevation in the shallow aquifer.
% When this water table intersects ground surface, surface runoff occurs.
% There is no limit to the infiltration capacity of the soil.
% Evapotranspiration maybe made depdendent on the extincktion depth.
% We may apply various tricks to influence the behavior of the boundary
% conditions such as more than one drainage level wiithin a single cell.
% The smallest conceivable model consists of only a single cell.
% It is sufficient to simuulate only half the strip of land between the
% ditch and its axis of symmetry at the center.
%
% TO 100905 101020
%
% Copyright 2009 2010 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

clear variables; close all;

basename='NPARK';  % Data for Noodererpark, near Utrecht Nettherlands

C_DEK_MIN=1;  % default minimum resistance cover layer 

%% Read meteo time series

fid=fopen('PE-92-00.txt','r');
A=fscanf(fid,'%d-%d-%d %f %f',[5,Inf])';
fclose(fid);

tne=[datenum(A(:,3),A(:,2),A(:,1)) A(:,[4 5])/1000,NaN(size(A(:,1)))]; % [t P N isSummer] % to mm/d

fprintf('Length of time series = %d\n',length(tne(:,1)));
%plot(tne(:,1),tne(:,2),'b',tne(:,1),tne(:,3),'r');

[YR,MONTH]=datevec(tne(:,1)');

tne(:,end)= MONTH>=4 & MONTH<=9;

clear A;

%% Get data from GGOR basis data spreadsheet (should be database)

data = dbfread(basename);

parnams = deblank({data.fieldname});

FID2 = data(strmatchi('FID2'    ,parnams,'exact')).values;
GP   = data(strmatchi('GP_ACT'  ,parnams,'exact')).values;  GP=round(100*GP)/100; % average phreatic head
ZP   = data(strmatchi('ZP_ACT'  ,parnams  )).values;  ZP=round(100*ZP)/100; % summer ditch level
WP   = data(strmatchi('WP_ACT'  ,parnams  )).values;  WP=round(100*WP)/100; % winter dithc level
L    = data(strmatchi('L'       ,parnams,'exact'  )).values;   L(L<20)=20; L =round(L);          % width of parcels
XC   = data(strmatchi('XC'      ,parnams,'exact'  )).values;   XC=round(XC);         % parcel center x
YC   = data(strmatchi('YC'      ,parnams,'exact'  )).values;   YC=round(YC);         % parcel center y
K    = data(strmatchi('K'       ,parnams,'exact'  )).values;   K =round(10*K)/10;    % parcel hor k
MU   = data(strmatchi('MU'      ,parnams,'exact'  )).values;   MU=round(100*MU)/100; % parcel Sy
C    = data(strmatchi('C_TRIW5' ,parnams    )).values; C(C>1000)=1000; C =round(C);          % parcel resistance (k_vert/d)
q    = data(strmatchi('QKW1'    ,parnams    )).values;   q(isnan(q))=0;        % upward seepage
phi  = data(strmatchi('PHI2'    ,parnams    )).values;   phi(abs(phi)>10)=0; phi=round(phi*1000)/1000; % head in regional aquifer
AHN  = data(strmatchi('AHNMEDCOR',parnams   )).values; AHN(abs(AHN)>10)=0; AHN=round(AHN*1000)/1000; % ground elevations
AREA = data(strmatchi('AREA'    ,parnams    )).values;
DDEK = data(strmatchi('DIKTE_DEK',parnams  )).values; % thickness of cover layer

GLGDBF = data(strmatchi('GLGNEW',parnams    )).values; GLGDBF = round(100*GLGDBF)/100; % rounded from database
GHGDBF = data(strmatchi('GHGNEW',parnams    )).values; GHGDBF = round(100*GHGDBF)/100; % rounded from database
GVGDBF = data(strmatchi('GVGNEW',parnams    )).values; GVGDBF = round(100*GVGDBF)/100; % rounded from database

GLGDBF(GLGDBF<-5) = NaN;
GHGDBF(GHGDBF<-5) = NaN;
GVGDBF(GVGDBF<-5) = NaN;

%% Set the physical model parameters

nRec = length(GP);

P = repmat(struct('FID2',NaN),nRec,1);
for i = nRec:-1:1
    
    % coordinates and depth
    P(i).FID2 = FID2(i);        % [Nr] identify this Xsection in the dbase file
    P(i).b    = L(i)/2;        % [m] half parcel size between water divide and center of ditch
    P(i).area = AREA(i);       % [m2]area repersented by parcel

    P(i).x    = XC(i);          % [m] parcel center x
    P(i).y    = YC(i);          % [m] parcel center y
    P(i).AHN  = AHN(i);
    P(i).D1   = DDEK(i);        % [m] thickness of top layer
    P(i).D1   = 6;              % in huidige GGOR vast
    P(i).D2   = 30;             % [m] default D of regional aquifer

    P(i).h_mean   = GP(i);     % [NAP] average phreatic head
    P(i).h_summer = ZP(i);     % [NAP] summer phreatic head
    P(i).h_winter = WP(i);     % [NAP] winter phreatic heat
    P(i).phi      = phi(i);    % [NAP] head in regional aquifer
    
    P(i).z0   = P(i).h_mean;          % [NAP] parcel ground elevation (in case of phreatic water)
    P(i).z1   = P(i).h_mean-P(i).D1;  % [NAP] default bottom of top aquifer
    P(i).z2   = P(i).z1-P(i).D2;      % [NAP] default bottom of second, regional aquifer'
    
    % Desired highest lowest and spring groundwater level from database,
    % computed with GGOR tool by Waternet, used here for comparison with
    % outcomes of this model
    P(i).GHGDBF   = GHGDBF(i) + P(i).AHN; % [NAP] GGG from databse (model Ouboter)
    P(i).GLGDBF   = GLGDBF(i) + P(i).AHN; % [NAP] GLG from databse (model Ouboter)
    P(i).GVGDBF   = GVGDBF(i) + P(i).AHN; % [NAP] GVG from databse (model Ouboter)
    
    % Prescribed seepage (upward positive) between phreatic and and second
    % aquifer. This was computed with an steady-state regional model
    P(i).q = q(i);            % [m/d] vertical (upward) seepage
    
    % Hydraulic parameters, conductivities of the phreatic and second
    P(i).hk1 = K(i);          % [m/d] hor hudraulic conductivity of top layer
    P(i).hk2 = 25;            % [m/d] default k of regional aquifer
    
    P(i).c   = max(C_DEK_MIN,C(i));  % [ d ] vert resistance of top layer
        
    % Vertical conductances or anisotropies of both aquifers (depends on
    P(i).vk1      = 0.5 * P(i).D1/P(i).c; % [ - ] vert anisotropy first layer, requires layvka to be 1 in the LAY worksheet
    P(i).vk2      = P(i).hk2;       % [ - ] default vertical cond regional aquifer, require layvka=1 in the LAY worksheet
       
    ss    = 1e-5;  % default speific storage coefficient
  
    P(i).sy1 = MU(i);      % [ - ] specific yield of first layer
    P(i).S1  = MU(i);      % [ - ] confined storage coefficient (MF2005 STORAGECOEFFICIENT option)
    P(i).ss1 = ss;      % [1/m] specific storage coefficient
    P(i).sy2 = MU(i);      % [ - ] specific yield of first layer
    P(i).S2  = ss*P(i).D2; % [ - ] elastic storage coefficient layer 2 (MF2005 STORAGECOEFFICIENT option)
    P(i).ss2 = ss;      % [1/m] default elastic storativity, all layers
    
    % Drains on ground surface to simulate surface runoff
    P(i).cdr   = 0.1;          % [ d ] default drain resistance
    P(i).zdr   = P(i).AHN;     % [NAP] default drain elevation    
    
    %% Ditch properties
    P(i).cdb   = 0.01;    % [ d ] default sludge resistance
    P(i).dw    =    2;    % [ m ] default ditch width    
    P(i).dd    =  0.6;    % [ m ] default ditch depth
    P(i).zdbot = P(i).h_mean-P(i).dd; % [ m ] ditch bottom elevation
      
    if P(i).zdbot>P(i).z1       % ditch does not cut through cover layer
        P(i).Omega1 = min(P(i).dw/2+P(i).dd,P(i).h_mean-P(i).z1);
        P(i).Omega2 = 0;
        P(i).w1     = P(i).cdb*P(i).D1/P(i).Omega1;  %+2/(pi*sqrt(P.hk1*P.vk1))*log(P.D1/P.Omega1*sqrt(P.hk1/P.vk1));
        P(i).w2     = 1e8;  %+2/(pi*sqrt(P.hk2*P.vk2))*log(P.D2/P.Omega2*sqrt(P.hk2/P.vk2));
  
        P(i).vk_ditch = min(P(i).vk2, 0.5*(P(i).zdbot-P(i).z1)/P(i).D1 / P(i).c); %+ 2/(pi*sqrt(P.hk1*P.hk2))*log(2*P.D2/P.dw *sqrt(P.hk1/P.vk2));
        
    else % ditch cuts into second regional aquifer
        P(i).Omega1= P(i).h_mean - P(i).z1;
        P(i).Omega2= P(i).dw/2+(P(i).z1-P(i).zdbot);
        P(i).w1    = P(i).cdb;                % because Omega1 = total depth of first aquifer
        P(i).w2    = P(i).cdb*P(i).D2/P(i).Omega2 + ...
            2/(pi*sqrt(P(i).hk2*P(i).vk2))*log(P(i).D2/P(i).Omega2*sqrt(P(i).hk2/P(i).vk2));
        
        P(i).vk_ditch = P(i).vk2;
    end
    
end

%% Cleanup
MINAREA = 1e3;            % [m2] minimum parcel area taken into account

% Select only those parcels with a minimum area if 1 ha and a width > 10 m
P = P([P.area]>=MINAREA); % Only use parcels > 1 hectare
P = P([P.b]>10);

%P=P(184);  % 1:250);
%P=P(1:50);
%P=P(1170);

save P; % need this for model testing in testmdl.m

clear parnams parvals GP ZP WP L XC YC K MU C q phi AHN  % clean up unnecessary variables

%% This changes the read data into something else: a test set
%  this allows testing the code
%  comment out to bypass the test run

%mktestset;  % change data into test set
verification
return;

%% Model grid
dx = 2;  % cell width choice
wd = 0.01;

xGr = [wd/2 0:dx:round(max([P.b]))]; % [m] width of model (max length of any parcel)
yGr = 0:1:length(P);                % each x-section gets its own row (CHANI=0)
zGr = [0;-1;-2];                    % default zGr values, will be adjusted to actual parcel

Z=NaN(Ny,Nx,length(zGr));
for iy=1:length(P)
    Z(iy,:,1)=P(iy).z0;
    Z(iy,:,2)=P(iy).z1;
    Z(iy,:,3)=P(iy).z2;
end

gr = gridObj(xGr,yGr,Z);

%% Zoning for later recognition of locations of boundaries

iDRN    =1;    % zone of drains
iDITCH1 =2;   % zone of ditch in 1st aquifer
iDITCH2 =22;  % zone of ditch in 2nd aquifer
iWVP2   =30;    % zone of regional aquifer
iWEL    =iWVP2;% zone where to dose seepage

%% Build IBOUND use zone number to identify zones

IBOUND=gr.const(1);

% some of these lines could be done witout loop, but this keeping all
% inside the loop provides a more direct comparison between individual
% lines
for i=1:length(P),
    IBOUND(i,xm>wd,1)=iDRN;    % location of drain
    IBOUND(i,xm<wd,1)=iDITCH1; % location of ditch in first  aquifer
    IBOUND(i,xm<wd,2)=iDITCH2; % location of ditch in second aquifer
    IBOUND(i,xm>wd,2)=iWEL;    % seepage will be divided over the second aquifer
    IBOUND(i,xm>P(i).b+eps,:)=0;   % make cells outside section inactive  
end

%% Build model arrays
[LAYparnams,LAYparvals]=getLayers(basename,Z);
LAYCON=LAYparvals(:,strmatchi('LAYCON',LAYparnams));

HK     = ones(size(IBOUND));
VK     = ones(size(IBOUND)); % use vertical conductivity not anisotropy
SY     = ones(size(IBOUND));
SS     = ones(size(IBOUND));
STRTHD = ones(size(IBOUND));

HK(:,:,1) =[P.hk1]' * ones(size(xm));  % kh of top aquifer
HK(:,:,2) =[P.hk2]' * ones(size(xm));  % kh of bottom aquifer
VK(:,:,1) =[P.vk1]' * ones(size(xm)); % kz of top aquifer
VK(:,:,2) =[P.vk2]' * ones(size(xm)); % kz of bottom aquifer

%% Elastic storage coefficients are used instead of specific storage
%  coefficients. Therefore STORAGECOEFFICIENT must be on and MF2005 must be
%  used (see worksheet MFLOW)
SY(:,:,1)   =[P.sy1]' * ones(size(xm));
SS(:,:,1)   =[P.S1 ]' * ones(size(xm));
SY(:,:,2)   =[P.sy2]' * ones(size(xm));
SS(:,:,2)   =[P.S2 ]' * ones(size(xm));

%% Make sure tha the ditch cell has a sufficiently highg conductivity
%  to prevent extra resistance within cell

for i=1:length(P)
     VK(i,xm<=wd,1)=P(i).vk_ditch;
end

%% mf2k has 1000 stress period limit, use mf2005 instead, see nam sheet

RECH = ones(1,1,length(tne(:,1)));
RECH(1,1,:)=tne(:,2)-tne(:,3);

%EVTR=ones(1,1,length(tne(:,1))); EVTR(1,1,:)=tne(:,3);

%% Start head

STRTHD(:,:,1)=[P.h_mean]' *ones(size(xm));  % start head of top aquifer
% STRTHD(:,:,2)=[P.phi]'*ones(size(xm));  % start head of bottom aquifer
STRTHD(:,:,2)=[P.h_mean]'*ones(size(xm));  % start head of bottom aquifer

%% Drains at ground surface to compute surface runoff

CDRN=NaN(size(IBOUND));              % allocate memory for drain conductance
CDRN(:,:,1)=(Dy*Dx)./([P.cdr]'*ones(size(Dx)));  % comptue drain conductance
CDRN(IBOUND==0)=NaN;                 % remove inactive cells
IDRN=find(~isnan(CDRN));             % get global incices of drain locations

hDRN=NaN(size(IBOUND));              % allocate memory
hDRN(:,:,1)=[P.z0]'*ones(size(Dx));  % [Ny*Nx] drain head is equal to z0, ground surface

%% Head boundaries in the left-most cell to simulate the ditch


%% GHB -- connection to ditch in first and second layer
CGHB1 = ([P.D1]./[P.w1])';     % comptatible with analytic solution
hGHB1 = [P.h_mean]';           % GHB heads in the left ditch only, whill change during stress periods
IGHB1 = find(IBOUND==iDITCH1); % global indices of ditch in 1st aquifer

CGHB2 = ([P.D2]./[P.w2])';     % comptatible with analytic solution
hGHB2 = [P.h_mean]';           % GHB heads in the left ditch only, whill change during stress periods
IGHB2 = find(IBOUND==iDITCH2); % global indices of ditch 2nd aquifer

%% Wells to simulate seepage from bottom aquifer or vice versa

% The specific upward positive seepage flux is given ("kwel").
% treat the seepage as uniformly distributed across the section.

% this is a bit of a quirk, take care qWEL matches IWEL
qwel= [P.q]'*ones(size(Dx));
IWEL = find(IBOUND(:,:,end)==iWEL); IWEL=IWEL(:);
qWEL = qwel(IWEL).*DX(IWEL).*DY(IWEL); qWEL=qWEL(:);

%% using the global incices compute the Layer Row Column indices necesary
%  to specify the boundary conditions in MODFLOE throuth the DRN, GHB and
%  WEL packages

LRCdrn  = cellIndices(IDRN ,size(IBOUND),'LRC'); % for surface runoff
LRCghb1 = cellIndices(IGHB1,size(IBOUND),'LRC'); % for ditches
LRCghb2 = cellIndices(IGHB2,size(IBOUND),'LRC'); % for ditches
LRCwel  = cellIndices(IWEL ,size(IBOUND),'LRC'); % for seepage

uDRN  = ones(size(IDRN )); % a column of ones for easyier multiplying below
uGHB1 = ones(size(IGHB1)); % same
uGHB2 = ones(size(IGHB2)); % same
uWEL  = ones(size(IWEL )); % same

%% specify the arrays for the mentioned boundary conditions so that
% mfLab can generate the input files for MODFLOW

[pernams,pervals]=getPeriods(basename);
NPER=size(pervals,1); % number of stress periods defined in worksheet PER

u=ones(size(MONTH));

iPer=1;

if tne(iPer,end)==1, hLR=[P.h_summer]'; else hLR=[P.h_winter]'; end

DRN  = [iPer*uDRN  LRCdrn  hDRN(IDRN) CDRN(IDRN) ];
GHB =  [iPer*uGHB1 LRCghb1 hLR CGHB1 ;...
        iPer*uGHB2 LRCghb2 hLR CGHB2 ];
WEL  = [iPer*uWEL  LRCwel  qWEL];
for iPer=2:NPER
    if tne(iPer,end)==1, hLR=[P.h_summer]'; else hLR=[P.h_winter]'; end
    GHB = [ GHB;  iPer*uGHB1 LRCghb1 hLR CGHB1 ; ...
                  iPer*uGHB2 LRCghb2 hLR CGHB2];     
    WEL  = [WEL;  -iPer ones(1,4)]; % same as in first period
    DRN  = [DRN;  -iPer ones(1,5)]; % same as in first period
end

save underneath P tne iDITCH1
