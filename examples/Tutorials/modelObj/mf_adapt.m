%% Tutorial of modelObj
%
% This file teaches and verifies the modelObj for use in manipulating models
% TO 120829

%% Introduction
% We will generate a Model consisting of an array of modelObj each containing
% a variable such as HK, VK, PEF, STRTHD, STCONC, IBOUND, RIV, GHB, grid,
% well etc, and will manipulate it by changing minimum layer thickness,
% joining and merging of layers and resampling the entire model to a
% completely new grid. We do this to demonstrate the various methods and to
% verify them where possible.

%% Initialize
clear variables; close all;

%% Set name of this model
basename='modelObjTutor';
save name basename

GREP = 'STRESS PERIOD';  % filters screen output of SEAWAT and MT3DMS

%% what is a modelOb and what ia Model in this context?
% a modelObj is a class whose instantionation (actual model) is as struct
% with three fields:
% * name  % the name of the variable like 'HK', 'gr', 'RIV' etc
% * var   % the actual variable, i.e.HK, gr, RIV etc
% * type  % the mflab type of the variable like '3Dlay', 'gridObj', 'stress'
%
% A Model is then just an array of modelObj instantiations.
%
% The class has members, i.e. a set of functions shared by all modelObj's
% and it has properties that have different values between the individual
% modelObj's.

%%
% show what is in class modelObj
modelObj

%% How to create modelObj's?
% See    help modelObj
% Call the modelObj constructor
%%
% Model(10) = modelObj(); %yields an Model consisting of 10 empty modelObj
%%
% in general call with arguments and limits:
% Model = modelObj(data,xlim,ylim,Ilay);
%
% where data is
% * the directory holding the arrays constituting the model in matfiles
%   such that the name of the mat file is the name of the varariable it
%   contains (say HK in HK.mat), and that the contents of the matfile is
%   the variable itself with this name, so HK, plus a second variable with
%   the name and the word type attached to it, so here 'HKtype' having as
%   content a string with the mfLab type of the variable, i.e. here '3Dlay'
%   meaning a 3D array with a layer for each model layer. Other types are
%   '3Dcbd' (3D array with a later for each confining bed), '3Dtime' for RECH, with
%   one layer for each stress period, 'stress' for RIV,GHB,DRN,CHD, 'zlist'
%   for variables like HANI with one value for each layer, 'gridObj' for
%   grid objects, 'wellObj' for an array of well objects and possibly more.
% * an array of modelObj, i.e. another model, which may have been saved
%   for instance by
%      save MyModel.mat Model
%   earlier and is now retrieved by
%      load MyModel
%   and passed to the modelObj.
%  Ix,Iy,Iz are the indices along the rows, columns and layers which are to
%   be extracted, in case a submodel is desired. These fields are optional,
%   if omitted the entire model is used.

%% Demo
% Generate a grid
xGr = -750:250:750; xm=0.5*(xGr(1:end-1)+xGr(2:end));
yGr = -400:100:400; ym=0.5*(yGr(1:end-1)+yGr(2:end));

z = [0   -25 -50 -75 -100 -125 -150 -175 -200]';
zm= [-20 -25 -30 -45  -50  -55  -60  -65 -100]';

Z = repmat(permute(interp1([xGr(1) mean(xGr) xGr(end)],[z zm z]',xm)',[3,2,1]),[length(ym),1,1]);

% These parameters are optional
MINDZ=0.1;
AXIAL=0;
LAYCBD=0;
%%
% the grid
gr = gridObj(xGr,yGr,Z,LAYCBD,MINDZ,AXIAL);
%%
% size of grid
gr.size

%% Generate HK and VK
hk= (1:gr.Nlay)'; %hk(3)=25/1000;
vk= hk;
HK = gr.const(hk);
VK = gr.const(vk);


%% Generate a Model consisting of 3 empty model array objects
Model(3) = modelObj();

%% Fill this model
Model(1).name = 'HK';
Model(1).var  =  HK;
Model(1).type = '3Dlay';

Model(2).name ='VK';
Model(2).var  = VK;
Model(2).type = '3Dlay';

Model(3).name = 'gr';
Model(3).var  = gr;
Model(3).type = 'gridObj';

%%
% Now model is an array of 3 modelObj with some properties filled, see
Model(1)
Model(2)
Model(3)

%%
% to see the names of the variables
{Model.name} %#ok
%%
% and their mfLab types
{Model.type} %#ok
%%
% or use display
display({Model.type});

%%
% save this model
save Model Model;

whos
%%
% clear the workspace
clear;
whos
%%
% load the model
load Model
whos

gr = Model.grid();
%%
% to cutout a submodel
xlim = gr.xGr([2 end-1]); ylim=gr.yGr([2 end-1]);
subModel = modelObj(Model,xlim,ylim);
%%
% size of the new submodel
grSub = subModel.grid();
%%
% size of old grid
gr.size
%%
% size of new grid
grSub.size

%% verify what happend to the variables
fprintf('Model(1):\n');
Model(1)
%%
fprintf('subModel(1):\n');
subModel(1)

%%
clear subModel grSub;

%% using methods
% Note that most methods work on all objects of the modelObj array, so even
% though Model is just an array of modelObj, the modelObj methods are
% generally called like this
% Model= Model.descr('another string');
% To obtain a new Model array (overwriting the old) with changed
% properties. The methods (in this case descr) is applied to all objects in
% the object array Model.

%%
% compute transmissivity and HK at specified cell indices or every where if omited
Ix = 3; Iy= 2;
[kD,Hk] = Model.transm(Ix,Iy);
%%
% same for vertical hydraulic resistance and VK
[c ,Vk] = Model.resistance(Ix,Iy);
%%
% these methods have been combined in the script verify to display these
% variables for all layers together for a selected Ix and Iy
Model = Model.descr('Initialized: Model with HK VK');
verify % see verify by typing:    edit verify

%%
% append a string to Model.description proper
Model  = Model.descr('Loaded: Full model Model.mat from disk');

%%
% Show the model and use the last description line in the title
Model.showBox(Model(1).description{1});

%%
% if you like you may save the model so far
% save Model00 Model

%% Convert CBD to LAY
% If the model was essentially made for the BCF package it can be changed
% to the LPF package by
% Model = Model.CBF2LPF();
%%
% for transport models it is strongly advised to remove confining beds and
% to change them into regular model layers:
Model = Model.removeCBD();
Model=  Model.descr('Removed: confining beds');
Model.showBox(Model(1).description{end});

verify;
 
%% Further grid manipulation
%%
% to adapt the mininum thickness of layers:
MINDZ = 5; layers = [1 5 6];
% given layers at least 1 MINDZ thick
Model = Model.mindz(MINDZ,layers);
Model=  Model.descr(sprintf('Enforced: layers [1 5 6] set ot %f',MINDZ));
Model.showBox(Model(1).description{end});

verify;
%%
% update grid in workspace
gr = Model.grid();

%%
% Generate a STRTHD array
dhdx         = 1/500;  % head gradient in x-direction
strthd       = 0 + dhdx*(gr.xm - gr.xm(1));

STRTHD       = repmat(strthd,[gr.Ny,1,gr.Nz]);
Model(end+1) = modelObj(STRTHD,'3Dlay');
Model        = Model.descr('Added:STRTHD');

ZLIST = permute(1:gr.Nz,[1 3 2]);
CHANI = ones(size(ZLIST));
Model(end+1) = modelObj(CHANI,'zlist');
Model(end+1) = modelObj(ZLIST,'zlist');
Model = Model.descr('Added: CHANI');
Model = Model.descr('Added: ZLIST');

ZTIM1 = repmat(ones(gr.Ny,gr.Nx),[1 1  1]);
ZTIM2 = repmat(ones(gr.Ny,gr.Nx),[1 1 13]);
Model(end+1) = modelObj(ZTIM1,'3Dtime');
Model(end+1) = modelObj(ZTIM2,'3Dtime');

%%
% Generate the IBOUND array
IBOUND       = gr.const(1);

%% Fix right-hand boundary of model using CHD (environmental head) and STCONC
iWest           = 2;  % arbitrary zone number not equal to 0 or 1
iEast           = 3;  % arbitrary zone number not equal to 0 or 1
CHDOPT          = 2;  % option to intepret specified heads as environmental heads in SEAWAT

% use IBOUND as zone array (optional)
IBOUND(:,end,:)    = iEast; % put zone number at right hand cells
IBOUND(:,  1,:)    = iWest; % put zone number at left  hand cells 
hE = STRTHD(IBOUND== iEast);
hW = STRTHD(IBOUND== iWest);

zoneVals = {iEast hE hE CHDOPT;  % one line per zone
            iWest hW hW CHDOPT};

Model(end+1) = modelObj(IBOUND,'3Dlay');
Model        = Model.descr('Added: IBOUND');

%%
% add to Model
ICBUND       = IBOUND;
Model(end+1) = modelObj(ICBUND,'3Dlay');
Model        = Model.descr('Added: ICBUND');

%% Test Arrays
ITEST = IBOUND; for i=1:size(ITEST,2), ITEST(:,i,:)=i; end
JTEST = IBOUND; for i=1:size(ITEST,1), JTEST(i,:,:)=i; end
KTEST = IBOUND; for i=1:size(ITEST,3), KTEST(:,:,i)=i; end
XTEST = IBOUND; for i=1:size(ITEST,1), XTEST(:,i,:)=i; end
YTEST = IBOUND; for i=1:size(ITEST,2), YTEST(i,:,:)=i; end
ZTEST = IBOUND; for i=1:size(ITEST,3), ZTEST(:,:,i)=i; end

Model(end+1) = modelObj(ITEST,'3Dlay');
Model(end+1) = modelObj(JTEST,'3Dlay');
Model(end+1) = modelObj(KTEST,'3Dlay');
Model(end+1) = modelObj(XTEST,'3Dlay');
Model(end+1) = modelObj(YTEST,'3Dlay');
Model(end+1) = modelObj(ZTEST,'3Dlay');


%%
% Generate CHD list
basename = 'modelObjTutor';
CHD = bcnZone(basename,'CHD',ITEST,zoneVals);

%%
% Add to the Model
Model(end+1) = modelObj(CHD,'stress');
Model = Model.descr('Added: CHD stress list');

%% get wells
%
[~,HK] = Model.transm();
well=gr.well(basename,HK,'wells');

Model(end+1) = modelObj(well,'wellObj');
Model = Model.descr('Added: wells');

%%
% Join some layers of the original model, while keeping total transmisivity
% and vertical hydraulic resistance.
%
% joinArray tells which of the old layers are merged in each new layer:
joinArray=[1 2 3 4 5 6 7 8;   % Old Layer Nrs
           1 2 2 3 3 3 4 4] ;  % New Layer Nrs

Model = Model.merge('join',joinArray,'z'); % coarsen grid in z-direction
Model = Model.descr('Joined: layers [1 2 3] and [4 5] too thin or irrelevant for the density problem.');
verify;
%%
% split ome layers to get more vertical detail
splitArray = [1 2 3 4;  % old layer numbers
              3 5 2 3]; % each split in this number of new layers

display(splitArray);
%%
% split it
Model = Model.merge('split',splitArray,'z');
Model = Model.descr('Split: layers 1-->3,3-->5');
verify;

%% Grid resampling
%%
% completely resample the grid using planes to keep in gr.Z and the number
% of subdivisions between each pair of kept planes. A plane is a layer
% surface, a model has one plane more than the number of layers and
% confining beds combined. Plane 1 is omitted, the last plane is implied
% but obliged if used. The number of planes equals the number of
% subdivisions specified.

%%
% Adapt the grid by a complete resampling
gr = Model.grid();
planesToKeep = [ 2  4 9 gr.Nz+1];
subdivisions = [ 5 10 10     20];
Model= Model.changeLayers( planesToKeep, subdivisions);

Model = Model.descr('Changed: layers using planes to keep and subdivisions');

Model.showBox(Model(1).description{end});
verify;

%%
peff = 0.35;
PEFF=gr.const(peff);

Model(end+1) =modelObj(PEFF,'3Dlay');
Model = Model.descr('Added: PEFF');
Model.showBox(Model(1).description{end});

%% Join rows in y-direction
% All rows in y direction are joined to obtain a cross section model in the
% x-direction. The join array is used to tell which of the old layers will
% be in each of the new layers. To join them all, the join array has numbers
% 1 through gr.Ny in the first rows (all rows) and all 1 in the seconde,
% meaning that all the top rows will become the first row of the new model.
joinArray= [1:gr.Ny;               % row numbers current grid
            ones(1,gr.Ny)];        % row numbers new grid

Model = Model.merge('join',joinArray,'y'); % join in y-direction
Model = Model.descr('Joined: all rows (y direction) to get a cross section representative for a given width.');
Model.showBox(Model(1).description{end});
Ix = 1; Iy=1;
verify;

%% Show the wells in this box, rotate the box (by hand) to see them
% The wells are plotted in 3D and are visible in the model by rotating the
% model in the figure by hand (use 3D rotation tool in top of figure
% window).

well = Model(strmatchi('wellObj',{Model.type})).var;

well = well.plotXY(gca); % plot 3D in existing figure

%% Show the wells
% Turn the model to better show the wells (by hand using the 3D rotation
% tool in the tool bar at the top of the figure screen

%%
% Get start concentration distribution using function
cFresh = 0;
cSea   = 1;
SigmaX = 0;
SigmaZ = 50;
zc     = -120;
xc     = 0;
STCONC = getInitialSalinity(gr,cFresh,cSea,xc,zc,SigmaX,SigmaZ);
%%
% Add STCONC to model arrays by appending a modelObj to Model.
% This way any new information can be added to the current Model, like even
% its computed heads and concentrations.

Model(end+1) = modelObj(STCONC,'3Dlay');   % and the variable type
Model = Model.descr('Added: Initial salinity distribution (STCONC).');

%% Set min thickness for somelayers

MINDZ = 3; layers = [4 8 12];
Model = Model.mindz(MINDZ,layers);
Model = Model.descr(sprintf('Changed: Minimum layer thickness for layers [4 8 12] to %g m',MINDZ));

Model.showBox(Model(1).description{end});
verify;

done;


%% Set min thickness for somelayers
% a more extreme change of  minimum layer thickness.
% Note that underlying layers will be compressed and end up having
% minimum layer thickness equal to gr.MINDZ.
% This small correction may cause the total transmissivity and resistance
% to change a small amount.
MINDZ = 15; layers = [25 30 40];
Model = Model.mindz(MINDZ,layers);
Model = Model.descr(sprintf('Changed: Minimum layer thickness for layers [4 8 12] to %g m',MINDZ));

Model.showBox(Model(1).description{end});
verify;

done;

%% Start working in the workspace
% Put all model arrays into workspace
unpack;

%% Show all actions done to generate the Model
% Note that all saved models can be retrieved via load Modelxx where xx the
% number used when saved. By calling Model.descr a print will be made of
% all actions taken to genereate the particular model:
Model.descr;

