%% Saltwater Intrusion Package, example 5

% Exampe 5, shape of brakish zone below a square island
% Consider a square island of 525 by 525 metrers. The origin of he
% coordinate system is chosen at the center of the island. Below the island
% is a 50 m thick aquifer with a hydraulic conductivity of 10 m/d. There is
% an infiltration on the island of 2 mm/s (specified with the recharge
% package). In addition, there is an area of 100 by 100 meter where there
% is a combined net discharge of 80 m3/d (the dark gray area in Figure 9a).
% The island is surrounded by an ocean with an equivalent freshwater level
% of 0.05 m; the vertical resistance to outflow in the ocean is c=2 days
% (Vcont=0.5).
% The domain is discretized into 21 by 21 cells of 50 by 50 m. A timestep
% of 50 days is used. The tip and toe tracking parameters are not too
% important (there is only a top where the brackish zone intersects the top
% of the aqufier); top and toeslope are 0.05; delzeta and zetamin are
% computed with equation 16. The ocean is represented by a strip of 150 m
% (6 cells) wide around the island (the light gray area in Figure 9a).
% Groundwater below the island consists of freshwater (nu=0), brackish water
% (nu=0.0125), and seawater (nu-0.025). Initially the brackish zone below
% the island is 10 m thick and extends from 25 to 35 m below the surface;
% below the strip respresenting the ocean, the surfaces bounding the
% brackish zone are at the top of the aquifer. Example 5 is summarized in
% Appendix E.


basename='swiex5';

%% Specifics of this example
K          = 10;    % hydraulic conductivity
peff       = 0.2;   % effective porosity
hOcean     = 0.05;  % ocean head specified
cOcean     =   2;   % resistance ocean bottom
iOcean     =   2;   % zone number ocean floor
iDischarge =   3;   % zone number discharging square
R          = 350;   % half width of island

%% Specify mesh

xGr=-525:50:525;
yGr=-525:50:525;
zGr=[0 -50];

gr = gridObj(xGr,yGr,zGr);

IBOUND = gr.const(1);

IBOUND(gr.ym<-R | gr.ym>R,:,1) = iOcean;
IBOUND(:,gr.xm<-R | gr.xm>R,1) = iOcean;
IBOUND(gr.ym>50 & gr.ym<200, gr.xm>-200 & gr.xm<-50,1) = iDischarge;

% Data arrays
HK     = gr.const(K);
VK     = gr.const(K);
STRTHD = gr.const(hOcean);
SSZ    = gr.const(peff);  % effective porosity called SSZ in SWI
ISOURCE= gr.const(0);     % zone of extraction or injection
ISOURCE(IBOUND==iOcean)=-3;
RECH   = gr.const(0.002); RECH(IBOUND==iOcean)=0;

%% WEL
Q = -80;
WEL = bcnZone(basename,'WEL',IBOUND,{iDischarge, Q/4});

ibound = IBOUND(:,:,1);
GHB = bcnZone(basename,'GHB',IBOUND,{iOcean, hOcean, gr.AREA(ibound==iOcean)./cOcean});

%% ZETA planes

ZETA.term{1} = gr.const(-25);
ZETA.term{2} = gr.const(-35);

ZETA.term{1}(IBOUND==iOcean) = 0;
ZETA.term{2}(IBOUND==iOcean) = 0;
