%% Example of a simple transport model
% Version 1 the basic flow model

% First choose a basename and make sure you have an excel workbook with
% this basename
global basename

basename = 'TransportModel2';
save('name','basename');  % save it for later retrieval

% you can get such a workbook by copying 'Parameters.xls' from
% mfLab/examples and rename it to 'TransportModel1.xls'

%% The basis for a transport model is always a flow model
% Version one is, therefore the basic flow model.

% But we need one without a confining bed and more refined

%% Uniform layer constants
hk   = [10, 0.005, 25];
vk   = hk/5;

%% Layer thickness and layer elevations
% aquifer thickness  25 and 65 m respectively
% aquitard thickness 10 m
% top of the model is assumed at z=100 and the water table at about z=0;
% doing this, makes sure that the water table is free to fluctuate within
% the top layer without becoming confined. The actual elevation of the
% water table will then be the result of the boundary conditions.
zGr = cumsum([0,-25,-10,-65]); zGr(1)=100;

%% model size is 2.5 km by 2.5 km in size, we chose uniform cell size of 50x50 m

% there is a well in the lower left corner of the model, so that the left
% and the front face of the model are planes of symmetry.
% the back face and right face of the model are considered no-flow
% the bottom of the model is also no-flow.
xGr = 0:50:2500;
yGr = 2500:-50:0;

%% Get a gridObject
LAYCBD = 0;  % non confining beds

gr = gridObj(xGr,yGr,zGr,LAYCBD);

%% Generate full 3D model grids for conductivities
HK     = gr.const(hk);
VK     = gr.const(vk);
STRTHD = gr.const(25); % high enough to make sure all cells have water
IBOUND = gr.const(1);
PEFF   = gr.const(0.35); % for the model layers
STCONC = gr.const(1);
ICBUND = IBOUND;       % cells to be computed

%% We will put drains on ground surface to allow it to pond.
% Generate a nice ground surface elevation first.

yDr = 12.5 * (1+erf((gr.ym-mean(gr.ym))/(sum(gr.dy)/3))); % s-curve along y
Z   = gr.Z;  % extracts a full 3D grid from gridObj

% generate a new top elevation for Z according to this s-curve shape
Z(:,:,1) = yDr*ones(1,gr.Nx);

% regenerate the grid, now using the new 3D z-array for elevations
gr = gridObj(xGr,yGr,Z,LAYCBD);

%% ===================================================================

% Gnerate model objects to that we can manipulate (refine) the model
model(8) = modelObj();
model(1).name ='HK';      model(1).var=HK;     model(1).type='3Dlay';
model(2).name ='VK';      model(2).var=VK;     model(2).type='3Dlay';
model(3).name ='IBOUND';  model(3).var=IBOUND; model(3).type='3Dlay';
model(4).name ='ICBUND';  model(4).var=ICBUND; model(4).type='3Dlay';
model(5).name ='PEFF';    model(5).var=PEFF;   model(5).type='3Dlay';
model(6).name ='STRTHD';  model(6).var=STRTHD; model(6).type='3Dlay';
model(7).name ='STCONC';  model(7).var=STCONC; model(7).type='3Dlay';
model(8).name ='gr';      model(8).var=gr;     model(8).type='gridObj';

planesToKeep = [1 2 3 4 ];
subDivisions = [ 7 3 13 ];
model = model.changeLayers(planesToKeep,subDivisions,'z');

% Unpack the model into the workspace
for i=1:numel(model)
    eval([ model(i).name '= model(i).var;']);
end

% ===================================================================

%% show the top surface that we now got
figure; hold on;
shading('interp');
surf(gr.xm,gr.ym,Z(:,:,1),'edgecolor','none');
view(3);

%% Let's put drains all over the top surface
iTop = 5;  % choose a zone number value (arbitraty if >0 and <>1)
cdr  = 1;  % choose an areal drainage [d]

IBOUND(:,:,1) = iTop; % set op of IBOUND equal to zone number

% get top ground elevation and op area for areal drains for chosen zone
% number iTop (we only have one zone in this case)
zoneValues = {iTop gr.Z(IBOUND==iTop) gr.AREA3(IBOUND==iTop)/cdr};

% Generate the list of drains specifying one drain in each top cell with
% the correct elevation and areal conductance value
DRN = bcnZone(basename,'DRN',IBOUND,zoneValues);

% In the PER worksheet
% choose only 25 years of stress periods,
% make sure transient=0 for steady-state computations.

% save the workbook by pressing <Ctrl-S>
% check if mf_adapt.m (this file) runs

% Then you can run this model by pressing
%mf_setup

% after the model finished successfully run
%mf_analyze

figure;   hold on;   view(3);   xlabel('x [m]');   ylabel('y [m]');   zlabel('z [m]');
gr.plotMesh('faceAlpha',0.15); % thin grey lines

save underneath

