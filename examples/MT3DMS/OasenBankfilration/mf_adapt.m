% Cross section showing functioning Large-Scale storage of desalinated
% water in Abu Dhabi. Photos are from field excursion in October 2010
% with the ISMAR conference.
%
% TO 071024 071127 110316 100327 110413
%
clear variables;
close all;

AXIAL=0;

basename='Oasen';

peff=0.35;
ss  =1e-5;
sy = 0.1;
cFresh=50/20000;  % relative salinity concentration
rivBotElev=-5; % [m NAP] elevation of river bottom
cRiver  = 30;  % [d] resistance of river bottom
cPolder = 50;  % [d] resistance of cover layer (eeklaag polder)
hPolder = -1.5; % maintained water level in the polder

%% Get configuration from excel workbook

xGr=-1500:10:1500;           % x grid lines [m]
yGr=[-0.5 0.5];              % y grid lines [m] make it a cross section of 1 m width;
zGr=[-1.5 -13 -23 -50   0:-5:-50 -50:-10:-100];   % add sufficient lines to make a sufficiently fine grid

gr = gridObj(xGr,yGr,zGr);

 HK       =mf_zone(basename,gr.xGr,gr.yGr,gr.zGr,'Config','Material','kh'); % horizontal conductivities
 [VK,Conf]=mf_zone(basename,gr.xGr,gr.yGr,gr.zGr,'Config','Material','kv'); % vertical   conductivities


%% We will use IBOUND to mark the cells that will be GHB and RIV
%  GHB for polder cells, representing cover layer resistance
%  RIV for river  cells
iPolder= 2;
iRiver = 3;

%% We will use RIV to fix river stage (river levels)

xRivL=Conf.xL(strmatchi('River',Conf.names)); % x of left  side of river
xRivR=Conf.xR(strmatchi('River',Conf.names)); % x of right side of river

% see that we overwrite IBOUND==3 where the river is
%% fixed head boundaries. Are interpreted as local point heads in fdm2dens
IBOUND=gr.const(1);
IBOUND(:,[1 end],2:end)=-1;  % left and right hand side are fixed head boundaries

IBOUND(:,:,1)                  = iPolder;  % use 3 for GHB
IBOUND(:,gr.xm>xRivL & gr.xm<xRivR,1)= iRiver;   % use 2 to indicate river cells in first layer

PEFF = gr.const(peff);
SS  =  gr.const(ss);
SY  =  gr.const(sy);

%% Start heads and concentrations

STCONC=gr.const(cFresh);
STRTHD=gr.const(0);

%% Extractions. 
well = wellObj(basename,'wells',gr,HK,'PER');

%% boundary conditions

ghbCond = gr.DX(IBOUND==iPolder).*gr.DY(IBOUND==iPolder)/cPolder;
rivCond = gr.DX(IBOUND==iRiver ).*gr.DY(IBOUND==iRiver) /cRiver;  % compute river bottom conductance

ghbVals = {iPolder,   hPolder,  ghbCond};
rivVals = {iRiver,  'RivStage', rivCond, rivBotElev};

GHB = bcnZone(basename,'GHB',IBOUND,ghbVals, cFresh);
RIV = bcnZone(basename,'RIV',IBOUND,rivVals,'concRiv');

ICBUND=IBOUND;

save underneath Conf well cFresh