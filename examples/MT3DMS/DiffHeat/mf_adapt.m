% MF_ADAPT -- MT3DMS benchmark. Simulation diffusion with and without
% advection using head gradient and well boundary conditions and two
% species, one of which is tempeature. Linear sorption is used.
% Dispersivitie are zero.
%
%  TO 100405 120505

basename='Diffusion';
GREP = 'STRESS';

% This example shows 4 cases of pure diffusion without groundwater flow.
%
% All four cases are easily verified with analytical solutions.
%
% Diffusion and conduction are computed using seawat (also possible
% directly in MT3DMS) without flow. They could alternatively be computed
% using the flow model so that head stands for concentration.
%
% To compare between mass diffusion and heat conduction, two species are
% used. To use two species, NCOMP and MCOMP maust be set to 2 in MT3D
% sheet.
% Also, to STCONC, STCONC{1} and STCONC{2} must be specified. Both cells
% contain a full 3D cell array with initial values for species 1 and
% species 2 respectively.
% If any point sources are used (SSM package), the point source
% concentrations for both species must be given (See documentation of SSM
% package in the MT3DMS manual, page 122.
%
% First stress period may be different from the rest to allow a block response
% by setting given concentration during first stress period only.
%
% Options WEL, CCC  (const conc cell) and MLS (mass loading source) are
% used also, see MT3DMS manual, SSM package, page 122

%% Parameters

% Parameter values are stored in the sheet DATA of the workbook basename
% (see actual basename above). They are imported here using getExceldata,
% specifying the workbookname, the sheetname and the vertical or horziontal
% order of the data labels. Vertical in this case.
[dnams,dvals]=getExcelData(basename,'data','vertical');

% Then we read out the required variables
L     = dvals(strmatchi('Lmdl',dnams),1);
dxmin = dvals(strmatchi('dxmin',dnams),1);
peff  = dvals(strmatchi('peff',dnams),1);   % effective porosity
rhos  = dvals(strmatchi('rhos',dnams),1);   % soilds density
rhow  = dvals(strmatchi('rhow',dnams),1);   % water density
rhob  = dvals(strmatchi('rhob',dnams),1);   % bulk density of porous medium
kh    = dvals(strmatchi('kh'  ,dnams),1);
kv    = dvals(strmatchi('kv'  ,dnams),1);
Cini  = dvals(strmatchi('Cini',dnams),1);   % initial conc
Tini  = dvals(strmatchi('Tini',dnams),1);   % initial temp
C0    = dvals(strmatchi('C0'  ,dnams),1);   % boundary conc
T0    = dvals(strmatchi('T0'  ,dnams),1);   % boundary temp
DmassD= dvals(strmatchi('DmassD',dnams),1); % diffusion coefficient
DtempD= dvals(strmatchi('DtempD',dnams),1); % diffusion coefficient
KdMass= dvals(strmatchi('KdMass',dnams),1); % distribution coefficient
KdTemp= dvals(strmatchi('KdTemp',dnams),1); % distribtuion coefficient
RetT  = dvals(strmatchi('RetT',dnams),1);   % retardation
RetC  = dvals(strmatchi('RetM',dnams),1);   % retardation
q     = dvals(strmatchi('q'   ,dnams),1);   % infiltration from left boundary if set

%% See if well package is on

[nm,txt] = xlsread(basename,'NAM');

wellIsOn = nm(strmatchi('WEL',txt(:,3)),3);

%% Mesh
xGr= 0:dxmin:L;
zGr= [-0.5 0.5];
yGr= [-0.5 0.5];

gr = gridObj (xGr,yGr,zGr);

%% Generate all other matrices
IBOUND = gr.const(99);
STCONC = {gr.const(Cini) gr.const(Tini)};
HK     = gr.const(kh);
VK     = gr.const(kv);
PEFF   = gr.const(peff);

%% Adaptation to boundaries

IBOUND(:,end,:)=-1; % fixed hand at right end of model

iL = 3;  % arbitrary but within IBOUND unique zone number of left boundary

%% Definition of the different cases (set thisCase to 1,2,3,4 or 5

thisCase= 6;

switch thisCase
    case 1 % temp and mass concentration fixed at x=0, no flow
        if wellIsOn, error('switch well package off'); end
        dhdx = 0; q=0;
        IBOUND(:,1,:)    = -iL;
        STCONC{1}(:,1,:) = C0;
        STCONC{2}(:,1,:) = T0;                
    case 2 % temmp and mass fixed at x=0, flow by given head gradient
        if wellIsOn, error('switch well package off'); end
        dhdx = -1/1000; q=0;
        IBOUND(:,1,:)    = -iL;
        STCONC{1}(abs(IBOUND)==iL) = C0;
        STCONC{2}(abs(IBOUND)==iL) = T0;                
    case 3,  % no gradient, using well instead of head gradient
        if ~wellIsOn, error('switch well package on'); end
        dhdx=0;
        IBOUND(:,1,:) = iL;
        zoneVal  = [ iL q ];  %each line[ZoneNr Value(=Q in caes of well))
        zoneConc = [C0 T0 ];  %each line conc for all species in this zone
        [WEL,PNTSRC] = bcnZone(basename,'WEL',IBOUND,zoneVal,zoneConc);
    case 4, % same as 3 but using time variant input (block pulse during first stress period)
            % see PER sheet columns C0 and T0
        if ~wellIsOn, error('switch well package on'); end
        dhdx=0;
        IBOUND(:,1,:) = iL;
        zoneVal  = [ iL q ];     % each line holds [ZoneNr Value(=Q in caes of well))
        zoneConc = {'C0' 'T0'};  % each line holds conc (column header in PER sheet)
        [WEL,PNTSRC] = bcnZone(basename,'WEL',IBOUND,zoneVal,zoneConc);
    case 5, % constant concentration cell during first stress peirod, under head gradient
        if wellIsOn, error('switch well package off'); end
        dhdx = -1/1000; q=0;
        IBOUND(:,1,:) = iL;
        zoneVal  = [ iL q ];     % each line holds [ZoneNr Value(=Q in caes of well))
        zoneConc = {'C0' 'T0'};  % each line holds conc (column header in PER sheet)
        [~,PNTSRC] = bcnZone(basename,'CCC',IBOUND,zoneVal,zoneConc);
    case 6, % mass loading at x=0 during first stress period and head gradient
        if wellIsOn, error('switch well package off'); end
        dhdx = -1/1000; q=0;
        IBOUND(:,1,:) = iL;
        zoneVal  = [ iL q ];     % each line holds [ZoneNr Value(=Q in caes of well))
        zoneConc = {'C0' 'T0'};  % each line holds conc (column header in PER sheet)
        [~,PNTSRC] = bcnZone(basename,'MLS',IBOUND,zoneVal,zoneConc);
    otherwise
        error('Set thisCase now<<%d>> to 1,2,3,4 or 5',thisCase);
end

STRTHD = dhdx*(gr.XMlay - gr.xm(end)); % guarantee flow to right

ICBUND = IBOUND;
%% Variable we may need in mf_analyze

save underneath thisCase Cini C0 Tini T0 peff DmassD DtempD KdMass KdTemp RetT RetC q dhdx kh

