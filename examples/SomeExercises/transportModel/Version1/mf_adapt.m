%% Example of a simple transport model
% Version 2 the basic flow model
% The same model as version 1, but with transport using MT3DMS

% First choose a basename and make sure you have an excel workbook with
% this basename
global basename

basename = 'TransportModel1';
save('name','basename');  % save it for later retrieval

%% Uniform layer constants
hk   = [10, 25];
vk   = hk/5;
vkcb = 0.001; % vertical k of confining ned

%% Layer thickness and layer elevations
zGr = cumsum([0,-25,-10,-65]); zGr(1)=100;

%% Model size is 2.5 km by 2.5 km in size, we chose uniform cell size of 50x50 m
xGr = 0:50:2500;
yGr = 2500:-50:0;

%% Get a gridObject
LAYCBD = 1;  % tell that there is a confining bed below layer 1

gr = gridObj(xGr,yGr,zGr,LAYCBD);

%% Generate full 3D model grids for conductivities
HK     = gr.const(hk);
VK     = gr.const(vk);
VKCB   = gr.constCB(vkcb);
STRTHD = gr.const(25); % high enough to make sure all cells have water
IBOUND = gr.const(1);


%% We will put drains on ground surface to allow it to pond.
% Generate a nice ground surface elevation first.

yDr = 12.5 * (1+erf((gr.ym-mean(gr.ym))/(sum(gr.dy)/3))); % s-curve along y

Z = gr.Z;  % extracts a full 3D grid from gridObj

% generate a new top elevation for Z according to this s-curve shape
Z(:,:,1) = yDr*ones(1,gr.Nx);

% regenerate the grid, now using the new 3D z-array for elevations
gr = gridObj(xGr,yGr,Z,LAYCBD);

% show the top surface that we now got
figure; hold on;
shading('interp');
surf(gr.xm,gr.ym,Z(:,:,1),'edgecolor','none');
view(3);

%% Put drains all over the top surface
iTop = 5;  % choose a zone number value (arbitraty if >0 and <>1)
cdr  = 1;  % choose an areal drainage [d]

IBOUND(:,:,1) = iTop; % set op of IBOUND equal to zone number

% get top ground elevation and op area for areal drains for chosen zone
% number iTop (we only have one zone in this case)
zoneValues = {iTop gr.Z(IBOUND==iTop) gr.AREA3(IBOUND==iTop)/cdr};

% Generate the list of drains specifying one drain in each top cell with
% the correct elevation and areal conductance value
DRN = bcnZone(basename,'DRN',IBOUND,zoneValues);

%% To do transport we need concentrations

% In the workbook sheet NAM
% choose MT3DMS
% further choose packages LST BAS DIS LPF, RCH, DRN OC PCG HDS BGT
% and for MT3DMS the packages FTL MLT6 + BTN ADV DSP SSM GCG
% (the rest is immaterial)

% In the PER worksheet
% choose only one stress period, length 365.24*25 (25 years)
% make sure transient=0 for steady-state computations.

% In the LAY worksheet
% LAYCON=1,  (free water table if head is below top of layer)
% CHANI=1,   (horizontal anisotropy for layer = 1)
% LAYVKA=1   (vertical LAYVKA is in fact VK, vertical K instead of vertical
%             anisotropy)

% save the workbook by pressing <Ctrl-S>
% check if mf_adapt.m (this file) runs

% Then you can run this model by pressing
%mf_setup

% after the model finished successfully run
%mf_analyze

%% ======= PARTICLE TRACKING =========

%% Do some particle tracking (see the exercies with mpath)
% for this to work make sure the model MODPATH is switched on in the NAM worksheet

% for particle tracking we need the effective porosity
PEFF  = gr.const(  0.35); % for the model layers
PORCB = gr.constCB(0.35); % for the confining beds

%% generate particles to start
% in our simple case we will just generate particles at the water table
% in some arbitrary column (since all columns are still equal)

zoneArray = gr.const(0);

% specify zone locations in zone array
zoneArray(:,10,1)=1;  % set zone 1
zoneArray(:,25,1)=2;  % set zone 2
zoneArray(:,40,1)=3;  % set zone 3

% we define 3 particle groups through the zones, which all pertain to zone
% array
zoneVals = {1,'name',basename,'placement',1,'IFace',6,'LineSpec','bo';
            2,'name',basename,'placement',1,'IFace',6,'LineSpec','ro';
            3,'name',basename,'placement',1,'IFace',6,'LineSpec','go'}; 

%% Particles
% Generate the mpath_partileGroupObj from which MODPATH can generate the
% starting locations of the particles.
% see help mpath_particleGroupObj for options

pGrp = mpath_particleGroupObj(gr,zoneArray,zoneVals); % particleGroups

%% Get the particles only to show them
pGrp   = pGrp.getParticles(gr);

%% Show particles in 3D
figure;   hold on;   view(3);   xlabel('x [m]');   ylabel('y [m]');   zlabel('z [m]');

gr.plotMesh('faceAlpha',0.15); % thin grey lines
pGrp.plot(); title('Particles starting points');

%% You can turn the graphic by hand to better view the particles

save underneath zoneVals

