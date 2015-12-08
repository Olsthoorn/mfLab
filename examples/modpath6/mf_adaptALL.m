% Example see USGS modpath Version 6 (2012), all examples named
% simulation1 through simulation5, see manual p36ff
% TO 130220

clear variables; close all

global basename

[~,basename] = fileparts(pwd);

%% This script is run by all examples, it defines the model grid and the model properties

%% Grid

%%
% All values are in feet as in the original example, but as long as
% a consisten system of units is applied this does not affect the answers.

xGr = 400*(0:25);  % xGrid coordinates, feet
yGr = 400*(0:25);  % yGrid coordinates, feet

%%
% The model has 2 aquifers, one phreatic 1 semi-confined (reprented by model cell layers)
% and 1 aquitards (confining beds). This implies 4 z-plane elevations
% have to be specified for this grid:
zGr = [400 270 220 200 100 0];

%%
% Define position of confining beds (below layer 1 and 2)
LAYCBD= 0;

%%
% Generate a grid object, which stores the grid and performs useful
% computations with the grid through its methods.
% (Type gr to see its contents)
gr=gridObj(xGr,yGr,zGr,LAYCBD);

%%
% Show the grid
%gr.plotGrid;

%% Aquitard and aquifer properties
% The layers are homogeneous in this example. Only one value is needed for
% each property and layer
kh   = [50 50 0.01 200 200];  % ft/d one value for each layer
kv   = [10 10 0.01  20  20];
por  = 0.3;
sy   = 0.2;
ss   = 0.0001;

%% Model arrays 3D

iRiv = 2;
%%
% Specify which cells have fixed heads?
IBOUND=gr.const(1);  % Any positive number indicatse normal cell
IBOUND(:,end,1)= iRiv;   % Fix head in left column layers 1 & 2
%%
% Initialize all starting head values at 0
STRTHD=gr.const(5000);

%%
% 3D property arrays
POR   = gr.const(por);
HK    = gr.const(kh);
VK    = gr.const(kv);
SY    = gr.const(sy);
SS    = gr.const(ss);

%% Get stress period info from worksheet PER in ex2.xls
[PERnams,PERvals] = getPeriods(basename);
NPER = size(PERvals,1);

%% Recharge

% See wokrsheet PER, becouase INRECH>0, RCH is read from RECH column,
% one value per stress period and constant across the model.
% recharge in first stress period as cell array

%% RIV
% Rivers may be generated in real-world coordinates of nodes of the river,
% with the arguments as required by the river package. 9i.e.
% stage,leakance,rivbottomElev.

%% RIV
% Then use the gridObj method gridObj/bcnZone to generate a list
% for the cells through which the drain passes, as required by MODFLOW:
% Each line of the list has [iPer iz iy ix head Leakance IFACE] Default iPer=1

RIV = bcnZone(basename,'RIV',IBOUND,{iRiv 320 1e5 315 6});

%% Wells
% Generate well objects using real-world coordinats given in sheet 'wells'
% of ex2.xls and the HK array to distribute the well extraction over the
% cells that are intersected by the well screen. The gridObj gr is used to
% place the wells in the actual grid. The extraction for each well is
% obtained from the column Q_?? or Q?? in the worksheet PER where ?? is the
% value of the number property of the well.
well = wellObj(basename,'wells',gr,HK,'PER');

%% MODPATH

%%
ZONE = gr.const(0);
ZONE(:,:,1)=1;      % zone in toplayer  = 1
ZONE(well.idx)=2;   % well  zone number = 2
ZONE(:,end,1)=3;    % river zone number = 3

%% The rest is done per simulation, it defines the particleStartingLocations
% and other modpath6 settings pertaining to the different simulations.