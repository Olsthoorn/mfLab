% MT3DMS Benchmark Problem see Zheng & Wang (1999), p140ff
% Two-dimensional Transport in a Radial Flow Field

%  TO 091114 091203 120417

basename='TwoD-Radial';
AXIAL=0;

% M3DMS Benchmark poblem: '2D radial flow field zee Zheng (1999) p140)'
% The test problem considered in this section concerns the two-dimensional
% transport of solute injected from a fully penetrating well. The problem
% is intended to test the accuracy of MT3DMS as applied to a radial flow
% system. The assumptions for this problem are:
% -- The injection rate of the well is constant
% -- The ambient groundwater velocity is negligible relative to the
%    velocity created by the injection
% -- The aquifer is homogeneous, isotropic, and infinite in areal extent.
% -- The flow fiel is steady state.

% Simulatiion duration 27d, see worksheet PER

dx=10;      % m
dy=10;    % m
dz= 1;      % m
peff=0.3;   % [] effective porosity, set in worsheet LAY, together with aT/aL = 1.0;
Q   =100;   % m3/d Injection rate
cInj=1;     % Injection water concentration
kh= 10;     % horizontal k
kv=10;      % vertical k, irrelevant if 1 layer

%% Mesh

NLAY=1; NROW=31; NCOL=31;

xGr = dx*(0:NCOL);
yGr = dy*(0:NROW);
zGr =-dz*(0:NLAY);

gr=gridObj(xGr,yGr,zGr);

%% Generate all other matrices
IBOUND=gr.const(1);  IBOUND(:,end,:)= -1; % fixed heads
ICBUND=gr.const(1);

%% Other inputs

HK       =gr.const(kh);
VK       =gr.const(kv);
PEFF     =gr.const(peff);
STRTHD   =gr.const(0);
STCONC{1}=gr.const(0);

iWel = 2; IBOUND(end,1,1)=iWel;

%% PNTSCR, to give wells their input concentration

[WEL,PNTSRC] = bcnZone(basename,'WEL',IBOUND,{iWel, Q},cInj);

save underneath cInj
