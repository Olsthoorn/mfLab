%% Script NHI (netherlands hydrologic instrument) to guide you through the
% NHI site (www.nhi.nu)
%
% TO 120430
%
% Trouble NHI data
% recharge.zip,c_laag.zip and startingheads.zip not connected to links
% S1 and S2 missing.
% for transient BCF HY is required
%
% dimensions of Z (cm ??)
% dimensions of river bottoms in cm?
% dimension of AHN in in m, however
% how to deal with the different surface waters? What is what?
% how deal with depth of surface water bottom?
% hat are sur wells? Dense well field in some parts of the Netherlands.
%
% The input file for the NHI modelA list of possible downloads is found under the page "Bibliotheek"
% Instead of downloading, you may want to copy the download locations to
% a spreadsheet first. (Point on one of them, right click and copy link,
% paste this link into spreadsheet).
%
% You may then run getNHIHdata.m to download the zip files. The mfiles uses
% the references you stored in the spreadsheet. Notice that the function is
% primitive as it containss som hard-wired references.Change them as
% needed.
%
% You see that they are all on the downloads page
% http://www.nhi.nu/downoads/
%
% You may navigate to it and investigate by yourself what's on.
% Clicking a file will download it to you computer.
% This page contained files that were missing on the Bibliotheek page or
% generated the error message "not found".
% (c_laag.zip, startingheads.zip and recharge.zip).
%
% If we know exactly which zip files contain what data files, we could
% write a script or function thad downloads an appropriate zipfile from the
% NHI site and, unpacks it, after which we reed them into Matlab and
% generate a 3D model array in full size of it. We can efficiently store
% such af file in a mat file with a proper name, such as kD.mat or C.mat
% Alternatively we first download all the zip files and unpack them to be
% sure we have everything we need.
% After this we start generating our model arrays.
% The advantage is that this procedure is more manageble, and is not
% depending on the quirks of the internet connection. The disadvantages is
% the large amount of storage requried on our own hard disk.
%
% Just to get familiar with the files and their contents we download all
% the zip files into subdirectory NHI/downloads1. And we also downlaod the
% unpacked files separately in NHI/downloads2. This requires double
% storage, but seems the only way to fairly easily verify the contents of
% each of the zip files.
%
% The mfiles geteNHIdata does this automatically.
%
% May may put the names of all files into our spreadheet (more than 100 and
% link them to their zip file.
%
% see NHI.xls
%
% once you have that, the model may be constructed array by array

fprintf('Setting some basic values for this model\n');
basename  = 'Amsterdam';
mfLabdir      = '../../../';
NHIdownloads  = 'examples/Geohydrology2/NHI/downloads2/';    % directory with NHI files
filesdir      = [mfLabdir NHIdownloads];
sheetNm   = 'files';
columnHdr = 'variable';

fprintf('basename   = %s\n',basename);
fprintf('filesdir   = %s\n',filesdir);
fprintf('sheetnm    = %s\n',sheetNm);
fprintf('columnHder = %s\n',columnHdr);

Nlay = length(dir([filesdir 'KD*.ASC']));
fprintf('The model has %d layer,\n',Nlay);

LAYCBD = ones(1,Nlay);     % all layers (except last) have a confining bed below them
MINDZ = 0.1;               % minimum layer depth
fprintf('LAYCBD = '); fprintf(' %d',LAYCBD); fprintf('. Meaning all aquifers have an aquitard belowh them (except the last)\n');
fprintf('MINDZ  = %.0f m. Minimum layer thickness used in NHI.\n',MINDZ);

%% Get Meta data

fprintf('Get meta data, using KD1.ASC file for the purpose.\n');
fid = fopen([ filesdir 'KD1.ASC'],'r'); meta = getNHImeta(fid); fclose(fid);

%% Set grid for NHI

fprintf('Compute grid coordinates and cell center coordinates of NHI model.\n');
xGrNHI = meta.XLLCORNER+(0:+1:meta.NCOLS)*meta.CELLSIZE;
yGrNHI = meta.YLLCORNER+(meta.NROWS:-1:0)*meta.CELLSIZE;
xmNHI  = 0.5*(xGrNHI(1:end-1)+ xGrNHI(2:end));
ymNHI  = 0.5*(yGrNHI(1:end-1)+ yGrNHI(2:end));

%% Model grid centered around pumping test  20x20 km mesh around pumping test
fprintf('Choosing a location to center your model corresponding''.\n');

location = 'Monument op de Dam';

fprintf('Chosen location: %s.\n',location);

[E, N] = kmlpath(location);

%%
% Change from wgs84 to Dutch national system coordiantes of this pumping station or location
fprintf('Computing Dutch RD coordinates from longitude latitude.\n');
[xw,yw] = wgs2rd(E,N);  % only the first ptest with this name
%%
% Round center coordinates to center of cell in NHI grid
fprintf('Rounding location to cell center.\n');
xC = 250*floor(xw/250)+125;
yC = 250*floor(yw/250)+125;

%% Define grid
% The grid is defined around the cell with xC,yC as its center and
% having ncells to the west, east, north and south and the same number of
% layers as the NHI model. This grid, therefore has 2*ncells in x-direction
% and also 2*ncells in y direction
ncells = 40; % 10 km

fprintf('Computing centers of %.0f by %.0f km grid around chosen location.\n',...
    2*meta.CELLSIZE*ncells*[1 1]/1000);
xLim = xC + (ncells+0.5)*[-meta.CELLSIZE meta.CELLSIZE];
yLim = yC + (ncells+0.5)*[-meta.CELLSIZE meta.CELLSIZE];

xGr = xGrNHI(xGrNHI>xLim(1) & xGrNHI<xLim(2));
yGr = yGrNHI(yGrNHI>yLim(1) & yGrNHI<yLim(2));

%% Get NHI indices of submmodel (i.e. indices of submodel within NHI grid
fprintf('Setting NHI indices of submodel to retrieve sublayer data.\n');
Ix = find(xmNHI>    xGr( 1      ) & xmNHI<    xGr(   end)  ); % xGr runs from low to high
Iy = find(ymNHI>min(yGr([1 end])) & ymNHI<max(yGr([1 end]))); % yGr runs from high to low

%% Using only these indices we extract the data arrays from the NHI arrays

fprintf('Getting transmissivities of all layers.\n');
TRAN   = makeArray(basename,sheetNm,filesdir,columnHdr,'TRAN',Ix,Iy);

save TRAN TRAN
%%
fprintf('Getting the elevation of the tops and bottoms of all aquifers\n');
fprintf('and confining bed in sequence.\n');
Z = makeArray(basename,sheetNm,filesdir,columnHdr,'Z',Ix,Iy);

Z(isnan(Z))= -0.40;  % AHN is NaN where open water --> about NAP -0.40 msl

%% Checking integrety of Z
%%
% 1: All Z except the top are in cm, so change Z(:,:,2:end) to m
fprintf('Converting all Z(:,:,2:end) from cm to m.\n');
fprintf('Top of model (AHN) is already in m.\n');
Z(:,:,2:end)=Z(:,:,2:end)/100;
%%
% 2: Top layer is not everywhere above its bottom, so make sure this is the
% case for all layers, not just the top layer, just to be sure
fprintf('Making sure that bottoms are at or below tops of layers\n');
for i=2:size(Z,3)
    dz = max(diff(Z(:,:,[i i-1]),1,3),MINDZ);
    Z(:,:,i)=Z(:,:,i-1)-dz;  % make sure elevation of first layer>=its bottom
end
%%
% 3: Check that no layer has negative thickness, print it
fprintf('Veryfying for each layer that bottom is below or at its top.\n');
DZ= abs(diff(Z,1,3));
for iz=1:size(DZ,3);
    fprintf('Layer %2d, DZ varies between %12.3f and%12.3f\n',iz,min(min(DZ(:,:,iz))),max(max(DZ(:,:,iz))));
end
%%
% 4: save
save Z    Z xGr yGr Ix Iy xGrNHI yGrNHI LAYCBD

%% Generate a grid for the submodel, use MINDZ=0;

fprintf('Generating a grid object (gr) for the submodel grid.\n');
fprintf('We use MINDZ=0, as is the case for NHI.\n');
gr    = gridObj(xGr, yGr, Z, LAYCBD, MINDZ);

%% Show cross sections with only the aquitards colored

Dsec = 5000; % Distance between successive cross sections in m
ND   = round(Dsec/meta.CELLSIZE); % same in cells

fprintf('Plotting cross sections along x-axis every %g km.\n',Dsec);
IyXS= 1:ND:gr.Ny;

for i=1:length(IyXS);
    figure; hold on; grid on; grey=get(gcf,'color');
    
    I=find(~isnan(gr.ZTcbd(IyXS(i),:,1)) & ~isnan(gr.ZBcbd(IyXS(i),:,end)));
    xLim=gr.xm(I([1 end]));

    gr.plotXSec(IyXS(i),[],'hlines',grey);
    title(sprintf('XS NHI, RD-y =%.0f',gr.ym(IyXS(i))-125));
    xlabel('x [m]'); ylabel('z [m]');
    set(gca,'xlim',xLim);
end

fprintf('Plotting cross sections along the y-axis, very %g km\n',Dsec/1000);
IxYS= 1:ND:gr.Nx;
%%
for i=1:length(IxYS);
    figure; hold on; grid on;
    I=find(~isnan(gr.ZTcbd(:,IxYS(i),1)) & ~isnan(gr.ZBcbd(:,IxYS(i),end)));
    yLim=gr.ym(I([end 1]));

    %% Plotting along y axis. Note the options 'smooth' 'on' trigging non-staircase lines
    gr.plotYSec(IxYS(i),[],'hlines',grey,'smooth','on');
    title(sprintf('YS NHI, RD-x =%.0f',gr.xm(IxYS(i))-125));
    xlabel('y [m]'); ylabel('z [m]');
    set(gca,'xlim',yLim);
end

%% Get vertical resistances of all confining beds

fprintf('Getting vertical hydraulic resistance of confining beds.\n');
C      = makeArray(basename,sheetNm,filesdir,columnHdr,'C',Ix,Iy');

fprintf('Changing them to VONCT=1/C for use by MODLFOW''s BCF package.\n');
VCONT = 1./C;

save VCONT VCONT

%% Get the recharge
fprintf('Getting the recharge\n');
RECH   = getNHIASC([filesdir 'RCH.ASC'],Ix,Iy);

save RECH RECH

%% Get the starting heads
fprintf('Getting starting heads of grid\n');
STRTHD = makeArray(basename,sheetNm,filesdir,columnHdr,'STRTHD',Ix,Iy);

save STRTHD STRTHD

%% Get the wells (for as far they are within the grid)
fprintf('Getting wells, only those within the grid.\n');
fname = getNHIfileNm(basename,sheetNm,columnHdr,'WELL'); 

%% Get wells
% normal wells, surwells(=irrigation wells?), kD at normal well sites and grid
[well,surwell,kD,WEL] = getNHISCD(gr,[filesdir fname],Ix,Iy);

save WEL WEL well surwell kD

%% Get the boundary conditions: RIV, GHB and DRN
fprintf('Getting boundary conditions (RIV, GHB,DRN)\n');

%% RIV first
fprintf('Getting RIV for interaction with main, primary, secondary and tertiary surface water.\n');
RIV0 = getNHIBCN({'RIVHMY'  'RIVCM'  'RIVBMY' },basename,sheetNm,columnHdr,filesdir,Ix,Iy);
RIV1 = getNHIBCN({'RIVHP1Y' 'RIVCP1' 'RIVBP1Y'},basename,sheetNm,columnHdr,filesdir,Ix,Iy);
RIV2 = getNHIBCN({'RIVHS1Y' 'RIVCS1' 'RIVBS1Y'},basename,sheetNm,columnHdr,filesdir,Ix,Iy);
RIV3 = getNHIBCN({'RIVHT1Y' 'RIVCT1' 'RIVBT1Y'},basename,sheetNm,columnHdr,filesdir,Ix,Iy);

%% Put RIV into layer that agrees with bottom depth
RIV0=setRIVdepth(RIV0,Z,gr.LAYCBD);
RIV1=setRIVdepth(RIV1,Z,gr.LAYCBD);
RIV2=setRIVdepth(RIV2,Z,gr.LAYCBD);
RIV3=setRIVdepth(RIV3,Z,gr.LAYCBD);

save RIV RIV0 RIV1 RIV2 RIV3

%% Get GHB, which are to set infiltration resistances, use RIV heads

fprintf('Getting GHB to be used for infiltration resistances:\n');
GHB0 = getNHIBCN({'RIVHMY'  'GHBCM' },basename,sheetNm,columnHdr,filesdir,Ix,Iy);
GHB1 = getNHIBCN({'RIVHP1Y' 'GHBCP1'},basename,sheetNm,columnHdr,filesdir,Ix,Iy);
GHB2 = getNHIBCN({'RIVHS1Y' 'GHBCS1'},basename,sheetNm,columnHdr,filesdir,Ix,Iy);

save GHB GHB0 GHB1 GHB2

%% Get drains Drainage, used for tile drainage only
DRN  = getNHIBCN({'DRNH' 'DRNC'},basename,sheetNm,columnHdr,filesdir,Ix,Iy);

save DRN DRN
%%
%%
figure; hold on; xlabel('x'); ylabel('y'); zlabel('z');

h=showLayers(gr,'gr.ZBlay',well,location);

%% Done
fprintf('Data for submodel %s ready.\n',location);

save grid gr


