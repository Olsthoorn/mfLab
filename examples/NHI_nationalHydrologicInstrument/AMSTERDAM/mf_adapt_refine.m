%% Extract submodel from NHI (netherlands hydrologic instrument)
%

%% Model grid will be centered around a point obtained from Google Earth
fprintf('Choose center of model using a Google Earth pin (i.e., a kml file copied to this directory)\n');

%% In this example we use the pin located in Lexmond Lexmond.kml

d = dir('*.kml');

[Path basename Ext]= fileparts(d(1).name);

fprintf('Chosen basename: %s.\n',basename);

[E, N] = kmlpath([basename Ext]);

%%
% Change from wgs84 to Dutch national system coordiantes of this pumping station or basename
fprintf('Computing Dutch RD coordinates from longitude latitude.\n');
[xw,yw] = wgs2rd(E,N);  % only the first ptest with this name

%%
% Round center coordinates to center of cell in NHI grid
fprintf('Rounding basename to cell center.\n');
xC = 250*floor(xw/250)+125;
yC = 250*floor(yw/250)+125;

load ../NHIdata/gr;  grNHI = gr;  

Ixc = hit(gr.xGr,xC);
Iyc = hit(gr.yGr,yC);


%% Define grid of 2xnCells in both the directions:
ncells = 80;

Ix = Ixc-ncells:Ixc+ncells;
Iy = Iyc-ncells:Iyc+ncells;

%% Using only these indices we extract the data arrays from the NHI arrays
gr = gridObj(grNHI.xGr([Ix Ix(end)+1]),...
             grNHI.yGr([Iy Iy(end)+1]),...
             grNHI.Z  ( Iy,Ix,:),...
             grNHI.LAYCBD,...
             grNHI.MINDZ);

%% TRAN (transmissivities)

load ../NHIdata/TRAN; TRAN=TRAN(Iy,Ix,:);

HY = TRAN./gr.DZlay; % HY required for convertable layers (LAYCON==3)


%% Show cross sections with only the aquitards colored
if 0

    ND   = 20;  %distance between cross section plots in cells

    fprintf('Plotting cross sections along x-axis');
    IyXS= 1:ND:gr.Ny;

    for i=1:length(IyXS);
        figure; hold on; grid on; grey=get(gcf,'color');

        title(sprintf('XS NHI, RD-y =%.0f',gr.ym(IyXS(i))-125));
        
        xlabel('x [m]'); ylabel('z [m]');

        gr.plotXSec(IyXS(i),[],'hlines',grey);
    end

    fprintf('Plotting cross sections along the y-axis');
    IxYS= 1:ND:gr.Nx;

    for i=1:length(IxYS);
        figure; hold on; grid on;
        
        title(sprintf('YS NHI, RD-x =%.0f',gr.xm(IxYS(i))-125));

        xlabel('y [m]'); ylabel('z [m]');
        
        gr.plotYSec(IxYS(i),[],'hlines',grey);
    end
    
end

%% Get other data
fprintf('Getting vertical hydraulic resistance of confining beds.\n');
load ../NHIdata/C.mat;  C=C(Iy,Ix,:); VCONT = 1./C;

fprintf('Getting the recharge\n');
load ../NHIdata/RECH;   RECH=RECH(Iy,Ix);

fprintf('Getting starting heads of grid\n');
load ../NHIdata/STRTHD; STRTHD=STRTHD(Iy,Ix,:);

%% Get the wells (for as far they are within the grid)
fprintf('Getting wells, only those within the grid.\n');
load ../NHIdata/well;
load ../NHIdata/WEL;
load ../NHIdata/kD;

I = find([well.ix]>=Ix(1) & [well.ix]<=Ix(end) & [well.iy]>=Iy(1) &[well.iy]<=Iy(end));

well=well(I);
kD  =kD  (I);

for i=1:length(I),
    well(i).ix  = well(i).ix-Ix(1)+1;
    well(i).iy  = well(i).iy-Iy(1)+1;
    well(i).idx = cellIndex(  well(i).ix,well(i).iy,well(i).iLay,gr.size);
    well(i).LRC = cellIndices(well(i).idx,gr.size,'LRC');
    kD  (i).ix  = well(i).ix;
    kD  (i).iy  = well(i).iy;
    kD  (i).idx = well(i).idx;
    kD  (i).LRC = well(i).LRC;
end

WEL = cutoutBCN(WEL{1},Ix,Iy);

fprintf('Getting RIV for interaction with main, primary, secondary and tertiary surface water.\n');
load ../NHIdata/RIV0
load ../NHIdata/RIV1
load ../NHIdata/RIV2
load ../NHIdata/RIV3

RIV = RIV0;   % RIV=[RIV0;RIV1;RIV2;RIV3];

RIV = cutoutBCN(RIV,Ix,Iy);

fprintf('Getting GHB to be used for infiltration resistances:\n');
load ../NHIdata/GHB0
load ../NHIdata/GHB1
load ../NHIdata/GHB2

GHB = GHB0;  %GHB = [GHB0;GHB1;GHB2]; 

GHB = cutoutBCN(GHB,Ix,Iy);

load ../NHIdata/DRN

DRN = cutoutBCN(DRN,Ix,Iy);

%% Show STRTHD in 3D

h=showLayers(gr,'STRTHD',well,basename);


%% Done
fprintf('Data for submodel %s ready.\n',basename);


%% Set the layer to 2 for all RIV, GHB and DRN that have a layer thickness < 1 m
%  Necesary to have rewetting working for thin layer above the water table.

if 1   % set 1 to 0 to demonstrate that this action is necessary
    coverThreshold = 1.0; %m

    fprintf('Lowering the layer number if layer thickness < %g m\n',coverThreshold);

    % Get global indices for RIV stress cells
    Idx = cellIndex(RIV(:,[4 3 2]),gr.size);
    
    % See which of them are in a cover layer thinner than coverThreshold
    tooThin = gr.DZlay(Idx)<coverThreshold;
    
    % Move these stresses one layer down
    RIV(tooThin,2)=RIV(tooThin,2)+1;
    
    % Communicate how many have been moved down
    fprintf('%d  RIV points lowered\n',sum(tooThin)); 

    % Same for GHB
    Idx = cellIndex(GHB(:,[4 3 2]),gr.size);
    tooThin = gr.DZlay(Idx)<coverThreshold;
    GHB(tooThin,2)=GHB(tooThin,2)+1;
    fprintf('%d  GHB points lowered\n',sum(tooThin)); 

    % Same for DRN
    Idx = cellIndex(DRN(:,[4 3 2]),gr.size);
    tooThin = gr.DZlay(Idx)<coverThreshold;
    DRN(tooThin,2)=DRN(tooThin,2)+1;
    fprintf('%d  DRN points lowered\n',sum(tooThin)); 
end


IBOUND=gr.const(1);  IBOUND(isnan(STRTHD))=0;

%% refine using a new grid

% Any new refined grid
xGrNew = %refineGrid using wells
yGrNew = %refineGrid using wells

% generate a temperary new grid with for zGr gr.zgr (vector of layer-average z-values
grNew  = gridObj(xGr,yGr,gr.zgr,gr.LAYCBD,gr.MINDZ);

% now refine one array or list after the other
% this could be done by getting them through a load in a function and
% then restoring them into the same mat file

Z      = gr.refine(grNew,Z):
TRAN   = gr.refine(grNew,TRAN);
VCONT  = gr.refine(grNew,VCONT);
STRTHD = gr.refine(grNew,STRTHD);

WEL    = gr.refine(grNew,WEL);
RIV    = gr.refine(grNew,RIV);
GHB    = gr.refine(grNew.GHB);
DRN    = gr.refine(grNew.DRN);

well   = gr.refine(grNew,well);

grNew = grNew(grNew.xGr,grNew.yGr,Z,grNew.LAYCBD,grNew.MINDZ);


%% Logical table within refine for arg1 and arg1
% arg1     arg2
gr         array

save underneath well kD
