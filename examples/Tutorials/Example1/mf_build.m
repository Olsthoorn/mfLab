%% Tutorial mfLab -- simple multi layer model
% Example see USGS Modflow 2000 manual, Open-File Report 00-92 with some
% changes to apply it with more convenient units (m and d instead of ft and sec)
% TO 090806 091129 120522

%% For this version, we make use of the LPF (Layer Property File) package, not BCF
% The only things that change when using LPF instead of BCF are that
% TRAN is replaced by HK.
% That VCONT is replaced by VKCB, the vertical hydraulic conductivities of
% the confining beds.
% Packages that are on: WEL DRN RCH

clear variables; close all;

%% Set the basename for this problem
% A basename is always given, as this will be the name of all files
% generated, which will only differ in their extension.
% Note that the workbook with the parameters must also have the name
% <<Basename.xls>>:
basename='example1';

%% Define a grid. Note that Modflow is always 3D.
xGr=0:50:1000;
yGr=0:50:1000;

%%
% The model has 5 layers, so 6 layer planes, here taken constant
zGr=[0 -50  -60 -90 -110 -160];

%% Specify LAYCBD if there are any confining beds
% The model has two confining beds, one below layer 1 and one below layer 2. This
% is specified by the vector LAYCBD. When the vector is defined in the
% workspace as is done here, this overrules the defintion in the worksheet
% LAY of the workbook. (See examples/Tutorial/gridObj).

LAYCBD = [1 1 0];

%%
% As can be seen from the vector LAYCBD, the model has three aquifers, i.e.
% layers as contrasted with the two CBD's (confining beds)

%% Generate a grid object,
gr=gridObj(xGr,yGr,zGr,LAYCBD);

figure; hold on; xlabel('x [m]'); ylabel('z [m]');
title('XSection with confining beds (CBDs) in grey');
gr.plotXSec(1,[],'hlines','r','vlines',[0.8 0.8 0.8],'smooth','title','XSection');
axis tight
%% Cells with fixed heads area specified through the IBOUND array.
% The IBOUND array is useful to define zones and to define which cells are inactive (IBOUND==0),
% with which cells will be computed (IBOUND>0) and which cells have a fixed head (IBOUND<0).

%% Generate the 3D boundary array (zoneArray) IBOUND using method const of the gridObj
% IBOUND is used as a zone array. Use a convenient value, for instance 99, to easily discern
% user-defined zones from a background like 0 or 1.
IBOUND=gr.const(99);

%%
% Fix the heads in the upper two aquifers at the left side of the model:
IBOUND(:,1,1:2)=-1;

%%
% Choose starting heads. Here they are all made zero.
STRTHD=gr.const(0);

%%
% Hydraulic conductivity vector
kh =[5; 8; 10];  % m/d

HK   = gr.const(kh);

%% Vertical conductivity for the confining beds
vkcb=[0.1; 0.1]; % m/d

VKCB = gr.const(vkcb);

%% Vertical conductivity
VK=HK/2;

%% Effective porosity
peff = 0.35;

PEFF = gr.const(peff);  % we need porosity for particle tracking
%% Recharge 3 mm/d

RECH= 0.003*ones(gr.Ny,gr.Nx,1);

%% Boundary conditions, DRAINS etc:
% Drains are generated here in real world coordinates with a leakance value
% expressed per m length of drain.
% x,y,Elev,Leakance values are specified for concrete points and are then
% interpolated linaerly by mfLab to add to nodes of the model

%% Specify the drains
drn=[ 50 375   -10 10;...
     500 375   -10 10];  % [x y Elev Leakance] assuming layer 1

%% Use the mehtod bcnLine of the grid object to put this drain in the model.

% The grid object method gr.bcnLine(data,'type') generates a list
% required by Modflow, one line per drain cell.
DRN=gr.bcnLine(basename,'DRN',drn,drn(:,4:end));

%%
% This list has the following contents:
% [iPer iz iy ix head Conductance]

%%
% mfLab requires the first item of each tuple to be the stres period
% number. This is to allow sorting the lines so that they can be presented
% to the goundwater code in the direct in the right order to be executed in
% sequence.
% The same method can be used for other boundary conditions such as CHD and RIV.

%% Lastly we need to make sure the wells are put correctly in the grid.
% Pertinent data for an undefined array of wells has been put in the worksheet
% wells of workbook <<basename.xls>>, where basename is example1 in this case.
% Pertinent data of the wells consist of at least the following list:
% [nr x y z1 z2 rw ]
% where z1 and z2 are the elevations of the top and the bottom of the screen
% of the well and rw is the well radius.
%
% Well objects are generated from this list in the workbook by means of the
% mehod well of the gridObj. To be able to succeed it needs the name of the
% workbook (basename) and the worksheet with the pertinent well
% data,'wells' in this case. It further requires the horizonal conductivity
% array to allow distributing the well extraction over the layers penetrated
% by the screen. This distribution is weighed by the relative
% transmissivity of each screen cell.
% The wells will look for their dynamic flow data in the columsn of the
% worksheet PER where stress period parameters are defined.The headers of
% the columns must therefore be sufficiently large.
% Generate both the well objects and the file WEL necesary for MODFLOW.

%% Starting points for particle tracking using MODPATH
well= wellObj(basename,'wells',gr,HK,'PER');

%%
% Show the grid
figure; xlabel('x [m]'); ylabel('y [m]');
gr.plotGrid;
%%
% Choose some zone Numbers
iZ1 = 66;
iZ2 = 67;

%%
% Select a contour (done by clicking on the screen with [xpond, ypond] = ginput
xpond = [ 522   483   393   370   310   271   344   453   575   603   699   734   699 522];
ypond = [ 317   382   370   317   346   490   654   694   692   545   466   358   268 317];

%%
% Plot this contour on the grid
plot(xpond,ypond,'g');

%%
% Find cells that are in the contour (with their center in the contour)
Inpond = inpolygon(gr.Xm,gr.Ym,xpond,ypond);

%%
% Add zone iZ1 to top layer of zoneArray IBOUND
IB1=IBOUND(:,:,  1); IB1(Inpond)=iZ1; IBOUND(:,:,  1)=IB1;
%%
% Add zone iZ2 to bottom layer of zoneArray IBOUND
IB2=IBOUND(:,:,end); IB2(Inpond)=iZ2; IBOUND(:,:,end)=IB2;

%% Specify the number of points and iFace for these zones
zoneVals = {iZ1, 5, 6 ; ...   % second arg is n, third is iFace
            iZ2, 1, 0};     % same for second zone

LOCATIONS = gr.startLoc(IBOUND,zoneVals);

figure; hold on; view(3); xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');

gr.plotMesh;

xyz  = gr.relloc2model(LOCATIONS);
for i=1:length(LOCATIONS)
    plot3(xyz{i}(:,1),xyz{i}(:,2),xyz{i}(:,3),[mf_color(i),'o']);
end
%%
% Anything not specifically necessary for the MDOFLOW model and, therefore,
% not automatically saved, can be saved in underneath.mat to be retrieved in mf_analyze
% after the models have finished.

save underneath gr drn well xpond ypond xyz
%%
% Before pressing <ENTER> check that the parameters in the workbook
% are correct:
% *Check that the correct model and packages are selected in the NAM worksheet.
% *Check that the stress periods are correct
% *Check layer variables are correct
% *Check that the boundary conditions are correct.
% *Check whether the run or individal stress periods are transient or not
% (parameter "Transient" in the PER worksheet.
% *Check the solver settings and so on.

%% Generate the input for MODFLOW and run MODFLOW by typing mf_setup
%% Analyse your results by running mf_analyze
%  After the groundwater model finished, we can analyze the results by running
%  mf_analyze<enter>
%  Adapt mf_analyze to your requirements.
