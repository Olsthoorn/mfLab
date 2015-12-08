%% The generic MODFLOW model in mfLab

%% How to run and visualize this model ?
%
% run this model by typing

% mf_setup; mf_analyze

%%
% on the command line

%% Illustrative variants of this model
% a) run this model with DRN and WEL off in the NAM sheet in workbook generic_steady
% b) then run it again with DRN on
% c) then run it again with WEL and DRN on

%% Aim of this model
% This mf_adapt implements a generic model, i.e. a model box that runs with
% minimal input. The idea is to generate this model and then to adapt it to
% you specific needs. This can be done in subsequent steps. This way, this
% generic model can be used as a startup of different models.

%% About the generic model itself
% This model is the example that has been used in MODFLOW manuals ever since
% 1988. However, we will use meters instead of feet and days instead of
% seconds. Further, we will explicitly use the full 3D LPF package for the
% cell-to-cell conductance compuations. We will not use the BCF, as this
% 'quasi-3D' package may be considered more or less obsolete, in practice.
% However, some love it and it can be used.
% For the LPF package to be used, it must be on in worksheet NAM in the
% accompanying workbook (generic_steady.xls).

%%
% The problem is visualized at the cover of the MODFLOW 2000 and 2005
% manuals. It consists of a groundwater system of 3 aquifers with 2
% confining beds/units separating them. Therefore, we have 5 layers in
% total, which have thickness 350, 50, 50 and 100 m respectively, with the
% top of the system at +200 m. This top is kind of arbitrily high, so that
% the top layers leaves room for a free water table, of which the elevation
% is computed and unknown a priori.

%%
% The model is 7500 by 7500 m, with cell sizes of 5000 5000 m. Unrealistically
% large cells are used, but they can be chosen smaller by the user.
% The model is closed  on all layers and all sides, except at the left boundary
% in layers 1 and 2, where its heads are fixed at zero.

%%
% The horizontal conductivities of the layers are 100 20 and 20 m/d
% respectively. The vertical hydraulic conductivities are taken the same as
% the horizontal values. Notice that the LPF packages forces to choose
% vertical conductivities, whereas the BCF package does not.
% The vertical resistance of the two confining units are 500 and 1000 d
% respectively.
% Model cell properties are specified as one value per layer in the LAY
% worksheet. (HK, VK, SS, SY, PEFF ...). They are automatically generated
% for this model and put in Matlab's workspace by using the mf_lay2mdl()function.

%% Stress periods
% Stress period information, including well discharge and any other
% information that may differ between stress periods is given in worksheet
% 'PER'. In this case there is one steady-sate stress period, with recharge 0.002 m/d.

%% Cauchy boundaries (DRN, GHB, RIV)
% A drain is active at y=37500, x = from 7500 to 47500 m with elevation
% of the drain increasing from 0 to 100 m along this line. The specific entry
% conductance is set equal to 1 [m/d], i.e. 1 m3/d per m  length per m head difference
% between the model cell and the drain.
% Specifying the drain in coordinates instead of cell numbers makes
% its spcecification independent of the grid.
% The DRN package must be on in the NAM worksheet.

%% Wells in general
% There is a number of wells active in the 3 layers. These are specified in
% the workbook in sheet 'wells'. The only info used of this sheet is the nr,
% the x, y, and z1, z2 of the specified wells. z1 and z2 are the top and
% bottom of the well screen given in any order.

%% Generating the generic model

%%
% First choose a basename. All input and output files will have the same
% basename with different extensions, that denote their contents. The
% extensions are defined for the different packages in the NAM worksheet.
basename = 'generic_steady';

%%
% save it for later retrieval
save('name','basename');


%% Make sure you have a workbook with the name [basename '.xls']

%% Specify the horizontal layout of the model
%  This is done by giving the grid line coordinates x and y direction.
%  Model coordinates are always parallel to to the rows and columns.
xGr = 0:5000:75000;
yGr = 75000:-5000:0;

%% 3D model properties arrays
% The 3D model arrays are conveniently generated using the previousl horizontal
% layout and one value per layer given for the different variables
% in the workhseet LAY. This concern the variabls HK (horizontal k), VK
% (vertical k), PEFF (effective porosity), SS (specific storage
% coefficient), SY (specific yield). VKCB (vertical conductivity of
% confining beds) and PORCB (effective porosity of confining beds) are
% taken from the columns PEFF and PERCB respectively for layers that are
% deemed to be confining beds as indicated by a 1 in the column LAYCBD.

%%
% SS and SY are only generated if the model is transient (in any of the stress
% periods). This info is obtained from the column 'transient' in worksheet 'PER'.

%%
% VKCB and PORCB are only generated if the model contains confining units.
% This info is obtained from the column LAYCBD in the worksheet LAY.

%%
% Layer thickness is obtained from the column D in worksheet LAY. The
% elevation of the top layer is missing and, therefore, is specified
% further down in this mfile.

%%
% The number of layers is deduced form the number lf lines in worksheet
% LAY, where those with 0 in column LAYCBD will be a model layer, and those
% with 1 in column LAYCBD will be a confining bed. Notice that the number
% of layer in MODFLOW is the number of layers with cells. Confining beds
% are not counted as layers but as confining beds that belong to a specific
% layers.

%%
% Additional detailing can be done later or separately by defining the
% mentioned arrays yourself directly or by changing them further down in this
% mfile.

%%
% To let this work make sure you have the columns in worksheet LAY with the
% corresponding headers: 'HK, VK, 'PEFF', 'SS', 'SY', 'D' (thickness), 'LAYCBD'
% Notice that VKCB is obtained from column VK and PORCB from column PEFF.


%% Generate the 3D model arrays of this generic model
% See the list in the command window showing what arrays have been
% generated.
mf_lay2mdl(basename,xGr,yGr);

%% gridObj gr is also generated
% mf_lay2mdl also generated a gridObj with the name gr. The gridObj takes
% care of all handling of the grid that may be necesary during modeling.
% See under mfLab/examples/tutorial/gridObj what you can do with a gridObj.

%% Change default top of gr.zGr to +200
% The generated gridObj has the default zTop of the model equal to zero.
% We can reset this to the desired value of 200 m by regnerating the grid
% object using its zGr and adding 200 m to it.
% See help gridObj.gridObj to get more info on generating a gridObj
gr = gridObj(gr.xGr,gr.yGr,gr.zGr+200,gr.LAYCBD);

%% Recharge
% The recharge will be taken from the PER worksheet (as long as it is one
% value per stress period, othewise we have to adapt it).
% Nothing needs to be specified here, becaus if the RCH package is ON in the
% NAM worksheet mfLab will grab the rch from the RECH column in the worksheet PER.
% To specify RECH yourself, generate an array of size [Ny,Nx,Nt] with the
% appropriate values. Here, Nt is the number of stress periods.
% Likewise, you may use the EVT package.

%% Wells
% The well locations and depth will be obtained from worksheet 'wells'
% This worksheet must at least containt the fields nr x y z1 z2. Any other
% fields will be ignored. In fact, they will be stored in the field
% well(i).Userdata.
well = wellObj(basename,'wells',gr,HK,{'PER','Q'});

%%
% The well flows will be obtained from worksheet 'PER', columns Q1, Q2 etc.
% as is specified with th e{'PER','Q'} argument. Note if there is only one
% column equal to 'Q' without well numbers, than all wells will use that
% column in the PER worksheet to get their flows (injection positive,
% extraction negative).

%% Fixed-head boundaries
% Fixed-head boundaries are specified through the IBOUND array by setting the
% values for the corresponding fixed-head cells equal to a negative number.
% Here we specify the left column of layers 1 and 2 to be fixed-headed.
IBOUND(:,1,[1 2]) = -1;
%%
% Cells can be excluded from the computation by setting their IBOUND value
% to 0.

%% Stresses, in this case DRN
% Heads and flows can also be specified indirectly using one of the stress lists
% WEL, DRN, GHB, CHD, RIV, CHD ..., corresponding to the packages with this name in
% the NAM worksheet.
% Defining stresses of the type 'WEL, DRN, RIV CHD etc all work the same
% way. Notice that using wellObj instead of WEL as was done above is generally
% more convenient as it is independent of the grid.

%%
% Here we specify a drain as a line object. The line coordinates are given
% as a [n by 3] array where the columns are the x, y and z coordinates
% respectively. Hence a drain, well, RIV etc can be an arbitrary polyline.
line = [ 7500 37500   0 ;
        47500 37500 100];

%%
% Here we generate the MODFLOW input for a DRN using the line as specified
% before and the specific conductance cDrn [m/d]
cDrn = 1; % [m/d]
DRN = gr.bcnLine(basename,'DRN',line,cDrn);

%%
% DRN is a cell array with per stress period the cells that pertain to the
% drain with their elevation and specific drain conductance. Note that these
% values may change between stress periods.
% Such a list contains the information for one cell on each line as follows:
%   [SP L R C   elevation  conductance]
% where
%   SP is stress period number
%   L R C are layer row and column of the cell
%   elevation of the bcn in cell
%   resistance of the drain in d/m

%%
% The label 'DRN' can be replaced by that of one the other stress pacakges.
% For other stress packages to worlk they have to be switched on in the NAM
% sheet.

%%
% The elevation of the line at the center of the cells will be used for the
% drain elevation. The drain conductance of a cell is the specific drain entry
% conductance multiplied by the length of the drain in that particlar cell.

%%
% Drains can only discharge, the drain discharge of a cell equals
% Q = (hCell-elevationDrn) * L cDrn [m3/d] if hCell-elevationDrn > 0

%%
% cDrn in the call to bcnLine may be
%     -- a single value, to be used for all cells pertaining to the drain
%     -- a list/array of values, one for each cell pertaining to the drain.
%     -- a string, which refers to the header in worksheet PER. This causes
%     the cDrn to vary between stress periods as defined in that column. 
%     In this case all cells of the zone get the same
%     value, but this value may change from stress period to stress period
%     according to the column in the PER worksheet.
%     To generate different values for each cell and stress period a more
%     advance input loop is necessary.

%% Finishing
% This completes the required input.
% run mf_adapt by typing its name

% mf_adapt

%%
% in the command window of Matlab to see of there are any errors.
% Verify the workspace to see that the required arrays are ok.
% Inspect them, including the grid.
% If ok, run MODFLOW (making sure that of the models in worksheet NAM
% only MDOFLOW2000 or MODFLOW2005 is switched on.
% Type

% mf_setup

%%
% to generate the input for MODFLOW and run it.
% When done, type

% mf_analyze
%%
% to visualize the results.





%% Same as above without any comments:
% Below the same lines as above are shown without comments. This shows that
% the input is quite compact. In subsequent models we wil uses this compact
% input explaining only what's new.
 
% basename = 'genericModel';  % name of this model
% save('name','basename');    % save it for later retrieval

% xGr = 0:5000:75000;    % grid x coords
% yGr = 75000:-5000:0;   % grid y coords

% mf_lay2mdl(basename,xGr,yGr); % generate model arrays (from LAY)

% gr = gridObj(xGr,yGr,gr.zGr+200,gr.LAYCBD); % shift Z upward by 200 m

% well = wellObj(basename,'wells',gr,HK,{'PER','Q'}); % get wells

% IBOUND(:,1,end) = -1;  % set fixed heads

% line = [ 7500 37500   0 ;   % coordinates of drain
%         47500 37500 100];

% DRN = gr.bcnLine(basename,'DRN',line,0.001); % drain input


%% TO 130614