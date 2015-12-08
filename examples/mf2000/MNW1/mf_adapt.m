% MNW1 benchmark Halford and Hansen 2002 USGS open-file report 02-293

% BUG in MNW1 -- the program crashes by segmentation fault if only one 
% MNW is used. Always use at least two MNW !

%  TO 110808 110822 121118

basename='MNW1';

%% Model data
hy     = 60;        % 60 ft/d
tran2  = 15000;     % ft^2/d
strthd = 200;       % ft
vcont  = 2e-3;      % 1/d
sf1    = [0.05; 0.0001]; % primary storage coefficient (whole aquifer)
sf2    = [0.1 ;    0.1]; % secondary storage coefficient (always specific yield LAYCON =2|3)
%% Zones
iDRN   = 20;
iCHD   = 21;

%% Mesh using table data
xGr = ( 0:  1:14)*2500;
yGr = (21 :-1:0 )*2500;
zGr = [200 50 0 -100];

LAYCBD= [1 0];

gr = gridObj(xGr,yGr,zGr,LAYCBD);

%% ARRAYS
IBOUND = gr.const(1);
HY     = gr.const(hy);    HY    = HY  (:,:,1);
TRAN   = gr.const(tran2);
VCONT  = gr.const(2e-3);  VCONT = VCONT(:,:,1);
SF1    = gr.const(sf1);    % primary storage coefficient
SF2    = gr.const(sf2);
STRTHD = gr.const(strthd);

%% Boundary conditions
IBOUND(13,6:end,1) = iDRN;
IBOUND( :,  end,1) = iCHD;

STRTHD(IBOUND==iDRN) = round(interp1([6 13],[131 128],6:13));
STRTHD(IBOUND==iCHD) = round(interp1([1 21],[139 119],1:21));

zoneDRN = {iDRN,  STRTHD(IBOUND==iDRN), 10000};
zoneCHD = {iCHD,  STRTHD(IBOUND==iCHD), STRTHD(IBOUND==iCHD)};

DRN = bcnZone(basename,'DRN',IBOUND,zoneDRN);
CHD = bcnZone(basename,'CHD',IBOUND,zoneCHD);

%% Wells
%MNW  = MNW1Obj(basename,'wells',gr,TRAN,'PER');
well = MNW1Obj(basename,'wells',gr,TRAN,'PER');

%% Anything for mf_analyze goes here
save Underneath iDRN iCHD
