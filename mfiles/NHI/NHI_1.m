%NHI_1 script to extract the data from NHI ASCII files in zipfiles at www.NHI.nu
%
% These files together form the cell-based data of the national hydrologic
% instrument, a currently 7 aquifer - 6 aquitard nation-wide groundwater
% model with 1200 by 1300 cells of resolution 250 by to 250 m.
% This file reads the ASCII files and save them the MATLAB way to gain
% speed in later useage and save diskspace to the extent posible while
% allowing later extraction of any part of the national model
%
% See the file
% extractFromNHI.m to guide you through the
% NHI site (www.nhi.nu)
%
% TO 120430
%
% for transient BCF HY is required
%
% dimensions of Z are in cm except for the land elevation (AHN-250) which
% is in m. All elevations are in national datum NAP (about mean sea level).
%
% Recharge is given as a 1200 by 1300 array for each cell and represents
% some average value. Transient data are not provide on the NHI site.
%
% Unclear at this pont how to deal with the different surface waters? What
% is what exactly?
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
% generate a 3D model array in full size of it.
%
% We store the data in one matfile NHIdata.mat using matfile function.
%
%
% Just to get familiar with the files and their contents we download all
% the zip files into subdirectory NHI/downloads1. And we also downlaod the
% unpacked files separately in NHI/downloads2. This requires double
% storage, but seems the only way to fairly easily verify the contents of
% each of the zip files.
%
% The mfiles geteNHIdata does this automatically.
%
% The names of all files are over 110. The spreadsheet NHI.xls links them to
% their original zip file and to the parmeter in the MODFLOW model.
%%

fprintf('Setting some basic values for this model\n');
basename  = 'NHI';              % also used for the xls file with the file names
filesdir  = './downloads2/';    % directory with NHI files
sheetNm   = 'files';            % sheet in basename.xls with names of unzipped files
columnHdr = 'variable';         % column in sheet with variable name connected to file

fprintf('basename   = %s\n',basename);
fprintf('filesdir   = %s\n',filesdir);
fprintf('sheetnm    = %s\n',sheetNm);
fprintf('columnHder = %s\n',columnHdr);

%% Determine the number of NHI layers (= aquifers) from #of files with KD* name
Nlay = length(dir([filesdir 'KD*.ASC']));
fprintf('The model has %d layer,\n',Nlay);

LAYCBD = ones(1,Nlay);     % all layers (except last) have a confining bed below them
MINDZ = 0.1;               % minimum layer depth, overruled by NHI files

fprintf('LAYCBD = '); fprintf(' %d',LAYCBD); fprintf('. Meaning all aquifers have an aquitard belowh them (except the last)\n');
fprintf('MINDZ  = %.0f m. Minimum layer thickness used in NHI.\n',MINDZ);

%% Get meta data from KD1.ASC. This is the same for all ASCII files that represent layer data

fprintf('Get meta data, using KD1.ASC file for the purpose.\n');
fid = fopen([ filesdir 'KD1.ASC'],'r'); meta = getNHImeta(fid); fclose(fid);

%% Set grid coordinates for NHI based on meta data just obtained

fprintf('Compute grid coordinates and cell center coordinates of NHI model.\n');
xGrNHI = meta.XLLCORNER+(0:+1:meta.NCOLS)*meta.CELLSIZE;
yGrNHI = meta.YLLCORNER+(meta.NROWS:-1:0)*meta.CELLSIZE;
Ix= 1:meta.NCOLS;
Iy= 1:meta.NROWS;

%% Open at mafile object to store the arrays binary on disk

NHIdata= matfile('NHIdata');

%% Get transmissivities

matObj.TRAN   = makeArray(basename,sheetNm,filesdir,columnHdr,'TRAN',Ix,Iy);

%% Get layer top and bottom elevations of layers and aquitards in sequence: Z(Nx,Ny,Nlac+Ncbd+1)

Z = makeArray(basename,sheetNm,filesdir,columnHdr,'Z',Ix,Iy);

%% Checking integrety of Z
%
fprintf('Converting all Z(:,:,2:end) from cm to m.\n');
fprintf('Top of model (AHN) is already in m.\n');
Z(:,:,2:end)=Z(:,:,2:end)/100;

fprintf('Making sure that bottoms are at or below tops of layers\n');
for i=2:size(Z,3)
    dz = max(diff(Z(:,:,[i i-1]),1,3),MINDZ);
    Z(:,:,i)=Z(:,:,i-1)-dz;  % make sure elevation of first layer>=its bottom
end

fprintf('Veryfying for each layer that bottom is below or at its top.\n');
DZ= abs(diff(Z,1,3));
for iz=1:size(DZ,3);
    fprintf('Layer %2d, DZ varies between %12.3f and%12.3f\n',iz,min(min(DZ(:,:,iz))),max(max(DZ(:,:,iz))));
end
%%
% generate grid
gr= gridObj(xGrNHI,yGrNHI,Z,LAYCBD,MINDZ);

%% Save the data
NHIdata.gr= gr;
NHIdata.Z = Z;
NHIdata.xGr = xGrNHI;
NHIdata.yGr = yGrNHI;

NHIdata.C       = makeArray(basename,sheetNm,filesdir,columnHdr,'C',Ix,Iy');
NHIdata. RECH   = getNHIASC([filesdir 'RCH.ASC'],Ix,Iy);
NHIdata.STRTHD  = makeArray(basename,sheetNm,filesdir,columnHdr,'STRTHD',Ix,Iy);

%% Get the wells (for as far they are within the grid)
fprintf('Getting wells, only those within the grid.\n');
fname = getNHIfileNm(basename,sheetNm,columnHdr,'WELL'); 

% normal wells, surwells(=irrigation wells?), kD at normal well sites and grid
[well,surwell,kD,WEL] = getNHISCD(gr,[filesdir fname],Ix,Iy);

NHIdata.well   = well;
NHIdata.surwell= surwell;
NHIdata.kD     = kD;
NHIdata.WEL    = WEL;

%% Get the boundary conditions (stresses): RIV, GHB and DRN
fprintf('Getting boundary conditions (RIV, GHB,DRN)\n');

%% RIV
fprintf('Getting RIV for interaction with main, primary, secondary and tertiary surface water.\n');
RIV0 = getNHIBCN({'RIVHMY'  'RIVCM'  'RIVBMY' },basename,sheetNm,columnHdr,filesdir,Ix,Iy);
RIV1 = getNHIBCN({'RIVHP1Y' 'RIVCP1' 'RIVBP1Y'},basename,sheetNm,columnHdr,filesdir,Ix,Iy);
RIV2 = getNHIBCN({'RIVHS1Y' 'RIVCS1' 'RIVBS1Y'},basename,sheetNm,columnHdr,filesdir,Ix,Iy);
RIV3 = getNHIBCN({'RIVHT1Y' 'RIVCT1' 'RIVBT1Y'},basename,sheetNm,columnHdr,filesdir,Ix,Iy);

% Put RIV into layer that agrees with bottom depth
NHIdata.RIV0=setRIVdepth(RIV0,Z,gr.LAYCBD);
NHIdata.RIV1=setRIVdepth(RIV1,Z,gr.LAYCBD);
NHIdata.RIV2=setRIVdepth(RIV2,Z,gr.LAYCBD);
NHIdata.RIV3=setRIVdepth(RIV3,Z,gr.LAYCBD);

%% GHB
% Used to set infiltration resistances, use RIV for heads

fprintf('Getting GHB to be used for infiltration resistances:\n');
NHIdata.GHB0 = getNHIBCN({'RIVHMY'  'GHBCM' },basename,sheetNm,columnHdr,filesdir,Ix,Iy);
NHIdata.GHB1 = getNHIBCN({'RIVHP1Y' 'GHBCP1'},basename,sheetNm,columnHdr,filesdir,Ix,Iy);
NHIdata.GHB2 = getNHIBCN({'RIVHS1Y' 'GHBCS1'},basename,sheetNm,columnHdr,filesdir,Ix,Iy);

%% DRN
% Used for tile drainage only
NHIdata.DRN  = getNHIBCN({'DRNH' 'DRNC'},basename,sheetNm,columnHdr,filesdir,Ix,Iy);

%% Show STRTHD in 3D

% For any part of the model coutout this part using indces and put this in
% a new grid (see how this in done in extactFromNHI.m in the example directories

%may be too big to show on smaller PCs
%h=showLayers(gr,NHIdata.STRTHD,NHIdata.well,'Countrywide model');

%% Done
fprintf('All NHI data save in file NHIdata.mat for retrieval by submodels\n');


