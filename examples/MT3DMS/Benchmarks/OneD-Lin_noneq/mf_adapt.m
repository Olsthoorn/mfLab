% mt3dms benchmark problem:
% One-dimensional Transport with Linear Nonequilibrium Sorption
% See  Zheng & Wang (1999) p137ff
%  TO 091112 091203

basename='OneD-lin-noneq';

% A feature of MT3DMS is the ability to simulate solute transport unnder
% either chemical or physical nonequilibirum. In this section (p137ff), MT3DMS
% is used to solve one-dimensional transport subject to linar but non-equilibrium
% sorption. The test problem is identialto that in the previous ssection on
% non-equilibrium sorption, except for the sorption constants, which are given
% here as the distribution coefficient Kd=0.933 cm3/g, and the first order mass
% transfer coefficient beta, ranging from zero to 20 /s.
%
% The entry concentration changes at time zero from zero to a constant value
% and returns to zero after a fixed time.
%
% The list of parameters to sovle the problem is given here:

kh      = 10;     % irrelevant because flow is given but required by Modflow
kv      =  0;     % irrelevant because only one layer, but required by Modlfow
peff    = 0.37;   % []
vx      = 0.1;    % cm/s,  % true velocity = vx/peff=
q       = vx*peff;
rhob    = 1.587;  % g/cm3                --- rhob in LAY worksheet
Kd      = 0.933;  % cm3/g                --- SP1 in LAY worksheet
beta    = [0 2e-3 1e-2 20]; % /s            --- SP2 in LAY worksheet
Co      = 0.05;   % mg/L
t0      = 160;    % s -- pulse duration  --- PER in PER worksheet

% The values here are not directly used except for the conductivities. We
% use the values given here in the worksheet LAY where we specify the
% sorption properties. t0 is used in the worksheet PER to define the length
% of the first stress period. The concentration Co is used in the
% definition of the PNTSRC for the SSM module of MT3DMS. vx, the seepage
% velocity is used below to compute the inflow in the left most cell, after
% which this value, then a WEL, was used in the WEL worksheet for the
% WELL-flow.

% We are dealing with only one type of process with different variables
% that we specifcy in the layer worksheeet. Therefore we can simulate
% the four beta cases in a single model run using four layers that are not
% connected as long as we take the vertical conductivity equal to zero.
% ISOTHM property in this case is 4.

% The worksheet LAY has 2 layers one numbered 0 and one numbered 1. One is
% for the Freundlich isotherm and one for the Langmuis isotherm. To switch
% between them change the layer numbers to make to correct one active and
% in worksheet MT3D chagne the value of ISOTHM to 2 for Freundlicht and to
% 3 for Langmuir.

%% The computational mesh
dx=0.16; dy=1; dz=1;

NROW=1; NCOL=70; NLAY=1;            % MESH definition

xGr=dx*(0:NCOL);
yGr=dy*(0:NROW);
zGr=dz*(0:NLAY);

gr=gridObj(xGr,yGr,zGr);

%% Generate all other matrices
IBOUND=gr.const(1);  IBOUND(:,end,:)= -1;
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

% boundary concentraion in column C1_1 in worksheet PER
[WEL,PNTSRC]= bcnZone(basename,'WEL',IBOUND,[iWest q*dz*dy],{'C1_1'});

