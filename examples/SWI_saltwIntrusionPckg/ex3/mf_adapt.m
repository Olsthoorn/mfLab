%% Saltwater Intrusion Package, example 4

basename='swiex3';

% Example 3,
% This example considers interface flow in a two-aquifer system. The plane
% of flow is 4000 m long and 2 m wide; both aquifers are 20 m thick and are
% separated by a leaky layers of 1 m thickness. The origin of the x,z
% coordinate system is chosen at the lower left-hand corner of the plane.
% The hydraulic conductivity of the top aquifer is 2 m/d and of the bottom
% aquifer 4 m/d; the effective porosity of both aquifers is 0.2. Flow is
% semi-confined at all times. The first 600 m of the top of aquifer 1
% represents the ocean floor; the freshwater head at the ocean bottom is
% 50 m. The leaky layer between aqfuiers 1 and 2 has a vertical resistance
% of flow of c = 100 days (Vcont=0.01). The vertical resistance to outflow
% in the ocean is c=50 days (Vcont=0.02). There is a compehensiver flow of
% 0.015 m2/d towards the ocean (eqaully distributed over both aqufiers at
% the buondary).
% Initially the interface between the frsh and sawater is straight and is
% at the top of each aquifer at x=500 m and at the bottom of each aquifer
% at x=900 m; i.e. the initial slope is 0.05. The dimensinless density
% nu=0.025.
% Each aquifer is discretized into 200 cells of 20 m long and the timestep
% is 1 year. Teh last cell in the top aquifer has a source of 0.01 m3/d,
% and the last cell in the bottom aquifer a source of 0.02 m3/d. The ocan
% bottom is modeled with the general head boundary (GHB) package. The head
% in the GHB cells is specified to teh quivalent frsh water head at the
% bottom of the ocean (50 m). The Cond value is computed as the Vont value
% (0.02) of the ocean bottom by the area of cell (40). The water in the GHB
% cells representing the ocean is salt. Hence, when the GHB cell acts as a
% source, the infiltrating water is zone 2. Whne the GHB acts as a sink,
% water that flows out of the cells is of the same type as the water at the
% top of the aquifer. Hence, the ISOURCE parameter is set to -2. The top
% and toe tracking parameters are a TOESLOPE of 0.02, a TIPSLOPE of 0.04, a
% DELZETA of 0.06 and a ZETAMIN of 0.006.

%% Specifics of this example
K         = [2 4];% hydraulic conductivity
cAquitard = 50;  % resistane between first and second aquifer
peff      = 0.2; % effective porosity
hOcean    = 50;  % ocean head specified
cOcean    = 50;  % resistance ocean bottom
iOcean    =  3;  % zone number of ocean water

%% Specify mesh
xGr=0:20:4000;
yGr=[0 2];
zGr=[41 21 20 0]; dz = abs(diff(zGr));
LAYCBD = [1 0];

gr = gridObj(xGr,yGr,zGr,LAYCBD);

%% Matrices of computation grid
IBOUND  = gr.const(1);

ibound = IBOUND(:,:,1);  ibound(:,gr.xm<600)=iOcean;

STRTHD  = gr.const(hOcean);
SSZ     = gr.const(peff);
ISOURCE = gr.const(1); ISOURCE(1,ibound==iOcean) = -2; 
VKCB    = gr.const(dz(2)/cAquitard); VKCB=VKCB(:,:,1);
HK      = gr.const(K');
VK      = gr.const(K');

%% 
well = wellObj(basename,'wells',gr,HK,'PER'); % see sheet wells and Q1 in sheet PER for Q

%% GHB
IBOUND(:,:,1)=ibound;

GHB    = bcnZone(basename,'GHB',IBOUND,{iOcean hOcean gr.AREA(ibound==iOcean)/cOcean},0);

%% ZETA planes
ZETA = XS([interp1([gr.xGr(1) 500 900 gr.xGr(end)],XS(gr.zGr([1 1 2 2])),gr.xm); ...
           interp1([gr.xGr(1) 500 900 gr.xGr(end)],XS(gr.zGr([3 3 4 4])),gr.xm)]);
