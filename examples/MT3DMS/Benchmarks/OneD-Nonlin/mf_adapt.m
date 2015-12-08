% mt3dms benchmark problem:
% One-dimensional Transport with Nonlinear or Nonequilibrium Sorptions
% MT3DMS benchmark problem, see  Zheng & Wang (1999) p135ff
%  TO 091112 120417

basename='OneD-nonlin';

% Grove and Stollenwerk (1984) presented a computer code for modeling
% one-dimensional advective-dispersive transport with nonlinear equilibium
% controlled sorption and ion-exchange. MT3DMS is used here to solve the
% same problem involving Freundlich and Langmuir isotherms.
%
% The entry concentration changes at time zero from zero to a constant value
% and returns to zero after a fixed time.
%
% The list of parameters to sovle the problem is given here:

kh          = 10;    % doesn't matter because we specify the flux and problem is 1D
kv          = 10;    % values are irrelevant in this 1d problem
peff        = 0.37;  % []
vx          = 0.1;   % cm/s, true velocity = vx/peff=
rhob        = 1.587; % g/cm3                --- rhob in LAY worksheet
Kfreundlich = 0.3;   % (ug/g)(L/mf)^alfa  --- SP1 in LAY worksheet
q           = vx*peff; % cm/s
alfa        = 0.7;   % []                 --- SP2 in LAY worksheet
KLangmuir   = 100;   % L/mg               --- SP1 in LAY worksheet
Slangmuir   = 0.003; % ug/g               --- SP2 in LAY worksheet
C0          = 0.05;  % mg/L
t0          = 160;   % s -- pulse duration  --- PER in PER worksheet

% The values here are not directly used except for the conductivities. We
% use the values given here in the worksheet LAY where we specify the
% sorption properties. t0 is used in the worksheet PER to define the length
% of the first stress period. The concentration Co is used in the
% definition of the PNTSRC for the SSM module of MT3DMS. vx, the seepage
% velocity is used below to compute the inflow in the left most cell, after
% which this value, then a WEL, was used in the WEL worksheet for the
% WELL-flow.

% We cannot use different sorption processes in different layers in a single run,
% because the ISOTHM property in the RCT module is set for the entire model for all
% species; only the spcific process variable values are set per layer or
% even per cell. Hence all species in the model are either subject to
% Freundlich or Langmuir sorption and, therefore te simulate both processe,
% the model has to be run at as many times. However, we may combine
% sorption with an irreversible reaction such as decay.

% Hence, contrary to the previous example, there is no value in have more
% than a single layer in the model.

% We will use 2 stress periods one of 160 seconds with Cleft=C0,
% and one wih 1500 seconds with Cleft=0; See bcnZone(..) combined with
% worksheet PER to get the C0-time information for the specified zone.

%% Mesh

dx=0.16; dy=1; dz=1;

NCOL=100; NROW=1; NLAY=1;

xGr = dx*(0:NCOL);
yGr = dy*(0:NROW);
zGr =-dz*(0:NLAY);

gr=gridObj(xGr,yGr,zGr);

%% Generate all other matrices
IBOUND=gr.const(1);  IBOUND(:,end,:)= -1;  % fix right hand side
ICBUND=gr.const(1);
STRTHD=gr.const(0);

STCONC{1}=gr.const(0);

HK    =gr.const(kh);
VK    =gr.const(kv);
PEFF  =gr.const(peff);

%% Wells at the west boundaries
% Using the WEST and EAST list it is now straigthforward to generate the
% list of wells, i.e. cells with given flow input. These are the cells at
% the west end of the model, as we fix the head at the right end.

iWest=2;

IBOUND(:,1,:)=iWest;

[WEL,PNTSRC]= bcnZone(basename,'WEL',IBOUND,[iWest q*dz*dy],{'C1_1'});

save underneath C0
