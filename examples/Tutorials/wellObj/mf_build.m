% Seawat Benchmark example
%  TO 090101 091207

basename='wellObjTutorial';
%BACKGROUND=1;

%% Parameters
L                   = 500;     % m half width of model
D                   = 20;      % m height and top of model
Z0                  = 0;       % ground surface elevation above datum [m]
peff                = 0.35;    % [ ] effective porosity
kh                  = 30;      % m/d hydraulic condictivity
kv                  = 30;
cSea                = 20000;
cFresh              = 100;     % freshwater concentraton
cFreshRel           = cFresh/cSea;
cSeaRel             = 1;       % kg/m3 brine concentration
rhof                = 1000;    % kg/m3 freshwater density
rhos                = 1025;    % kg/m3 brine density
hCorner             = 150;     % m (point water heads)
HM                  = 5; grad = 0; % [m/m] head in center and head gradient with x dh/dx


%% Mesh
x                   = sinespace(0,L,50,pi/2);
xGr                 = [-x x];
zGr                 = [Z0:-0.5:Z0-D Z0-D];
yGr                 = [0 1];

gr=gridObj(xGr,yGr,zGr);

%% Generate all other matrices using gridObj/const method
IBOUND              = gr.const(99);
ICBUND              = gr.const(99);
STRTHD              = gr.const(0);
STCONC              = gr.const(cSea);
HK                  = gr.const(kh);
VK                  = gr.const(kv);
PEFF                = gr.const(peff);

%% Generate STRTHD given a groundwater gradient
STRTHD(:,:,1)       = HM+gr.xm*grad;

izoneCHD=2;       % Zone number for CHD

CHDDENSOPT=2;     % Option to set hydrostatic salt boundary. see Seawat manual

%% Get the number of stress periods
% This number allows adding the PNTSRCs from different boundary
% conditions, i.e. wells and constant heads CHD to be merged.
[~,~,NPER]=getPeriods(basename);

%% Generate well objects, the WEL array that is necessary for the WEL package and the PNTSRC array for the BTN package 
% The informmation originates from the workbook and the worksheet in which the
% pertinent well information is stored ('basename/WELLS' in this case).
% The conductivity array is used to
% estimate the fractions extracted from the cells that are penetrated by
% the same well screen. Well flow information comes from the PER worksheet,
% i.e. the columns Q_1, Q_2 etc.
well = wellObj(basename,'wells',gr,HK);

%% Use the IBOUND array as a zone array
IBOUND(:,[1 end],:)=izoneCHD;  % put zone into IBOUND array

%% Generate the arrays or cells necessary to feed MODFLOW and MT3DMS/SEAWAT

zoneVals = {izoneCHD, STRTHD(IBOUND==izoneCHD), STRTHD(IBOUND==izoneCHD), CHDDENSOPT}; % only one zone

[CHD,PNTSRC]=bcnZone(basename,'CHD',IBOUND,zoneVals,cSea);

%% Save anything that may be useful in mf_analyze
save underneath well rhos rhof L cFresh cSea cFreshRel cSeaRel gr peff kv kh
