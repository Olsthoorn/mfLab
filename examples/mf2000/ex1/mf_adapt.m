%% EX1 -- Example see USGS Modflow 2000 manual
% First example described in the USGS Open-File Report 00-92
% for modflow 2000.
% This file generates the required model arrays for this example.
% These arrays will be place in the Matlab workspace.
% Additional information will be obtained from workbook ex1.xls.
%
% Run mf_setup from the command line to generate modflow input files
% and run modflow.
%
% Visualize output by running mf_analyze thereafter, once modflow is ready.
%
% TO 090806 091129

clear variables; close all; % ensure memory is clean and there are no remaining figures

basename='ex1'; % name of this model, basename of all its files

%% Model grid

% specify grid coordinates
xGr    = 0:5000:75000;                     % [ft] xGrid coordinates
yGr    = 0:5000:75000;                     % [ft] yGrid coordinates
zGr    = [200 -150  -200 -300 -350 -450];  % [ft] Plane elevation vector

LAYCBD = [1 1 1]; % tell which layers have a confining bed below them

gr = gridObj(xGr,yGr,zGr,LAYCBD);          % generate grid object

%% Model arrays

IBOUND = gr.const(99);              % Boundary array, by >0 indicas cells to compute
IBOUND(:,1,1:2) = -1;               % Overwrite left column for layers 1 & 2, by <0 --> fixed heads

STRTHD = gr.const(0);               % [ ft  ]starting heads are all zeros

HY     = gr.const(1e-3);            % [ft/s ] kh of water-table layer

TRAN   = gr.const(1e-3*[1 1 2]');   % [ft2/s] transmissivity of layers 2 and 3

VCONT=gr.const(1e-8*[2 1 0]');      % [ 1/d ] Leakance of confining beds

RECH={3.0e-8*ones(gr.Ny,gr.Nx,1)};  % [ft/s ] recharge in first stress period

%% Boundary conditions

drn  = [5000 37500     0 1/5000;... %
        50000 37500  100 1/5000];   % [x y Elev Leakance] assuming drain in layer 1
 
DRN  = gr.bcnLine(basename,'DRN',drn,drn(:,end)); % DRN boundary array

well = wellObj(basename,'wells',gr,TRAN,'PER');   % geneate well objects

%% Saving non standard variables required in mf_analyse

save underneath drn %  anything else you might need in mf_analyze
