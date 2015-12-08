%% Saltwater Intrusion Package, example 4

basename='swiex4';

%% Upconing of an interface below a pumping well in a two-aquifer system
% Consider 3D flow in the same 2-aquifer system as in example 3.
% The area of interest extends 1300 m in the x-direction,
% and 1000 m in the y-direction. The upper aquifer is bounded on top by the
% ocean floor in the west (Figure 8z); a 300 m wide section of ocean floor
% is modeled. The flow entering the domain in the east is 0.12 m2/d, so that
% the total flow entering the system is 120 m3/d. A well of discharge Q=70
% m3/d is situated in the top aquifer at (x,y)=(975,525).
% Initially, the interface betwen the fresh and seawater is straight and is
% at the top of each aquifer at x=250 m and at the bottom of each aquifer
% at x=450 m; i.e. the initial slope is 0.1. The dimensionless density of
% the seawater is nu=0.025.
% Each aquifer is discretized in 20*26 square cells of 50 by 50 m; the well
% s at the center of cell (10,20). The timestep is one year. The top and
% toe tracking parameters are a TOESLOPE and TOPSLOPE of 0.05, a DELZETA of
% 0.25 and a ZETAMIN of 0.025. Example 4 is summarized in table 4; all
% input files are in Appendix D (of the SWI manual, part 1).
% Teh position of the interface after 50 years is to be shown as wel as the
% position after 400 years.

%% Specifics of this example
K         = [2 4];% hydraulic conductivity
cAquitard =  50;  % resistane between first and second aquifer
peff      = 0.2;  % effective porosity

hOcean    = 50;   % ocean head specified
cOcean    = 50;   % resistance ocean bottom
iOcean    = 3;    % zone number of ocean floor
iEastTop  = 4;    % zone number of east boundary top aquifer
iEastBot  = 5;    % zone number of east boundary lower aquifer

%% Grid
xGr=0:50:1300;
yGr=0:50:1000;
zGr=[41 21 20 0];
LAYCBD = [1 0];

gr = gridObj(xGr,yGr,zGr,LAYCBD);

%% Matrices of computation grid
IBOUND  = gr.const(1);
HK      = gr.const(K);
VK      = gr.const(K);
VKCB    = gr.DZ(:,:,2)/cAquitard;

STRTHD  = gr.const(hOcean);
SSZ     = gr.const(peff);
ISOURCE = gr.const(1);  ISOURCE(:,gr.xm<300,1) = -2;

%% wells and inflow from right implemented as wells
well = wellObj(basename,'wells',gr,HK,'PER');

%% use wells for east side of model to specify flow
IBOUND(:,end,1) = iEastTop;
IBOUND(:,end,2) = iEastBot;
Qtop = 2;
Qbot = 4;

zoneVals = {iEastTop, Qtop; iEastBot, Qbot};

WEL    = bcnZone(basename,'WEL',IBOUND,zoneVals);

%% GHB
ibound = IBOUND(:,:,1);
ibound(:,gr.xm<300) = iOcean;
IBOUND(:,:,1) = ibound;

GHB    = bcnZone(basename,'GHB',IBOUND,{iOcean, hOcean, gr.AREA(ibound==iOcean)/cOcean},0);

%% ZETA planes

ZETA=XS([interp1([gr.xGr(1) 250 450 gr.xGr(end)],[41 41 21 21],gr.xm); ...
         interp1([gr.xGr(1) 250 450 gr.xGr(end)],[20 20  0  0],gr.xm)]);

               
ZETA=repmat(ZETA,[gr.Ny,1,1]); % extend over all rows
