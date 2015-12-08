%NHI extrac NHI data from the ASCII files in mfLab/examples/NHI/NHIascii
%
% These files on their turn were obtained though unzipping the files in
%    mfLab/examples/NHI/NHIzipfiles. The zipfiles
% which were downloaded from the official NHI site (www.NHI.nu/bibliotheek).
% TO 120401

%% The NHI files
% These NHI files consitute the cell-based data of the national hydrologic
% instrument, a currently 7 aquifer and 6 aquitard nation-wide finite
% difference groundwater (and surface water) model. The exent of the model
% is 1200 by 1300 cells with constant cell resolution 250 by to 250 m.
%
% This script,
mfilename
%%
% reads the ASCII files and saves them as MATLAB mat files to gain maximum
% speed in later usage while using least diskspace.
% The total size of the ASCII files so reduces to around 250 Mb.
% To my knowledge this is the most efficient way of storage we can get.
%
% The stored files are generally either 3D arrays of some parameter or
% cell-based lists specifying stresses for the MODFLOW model. The stresses
% are for the WEL, GHB, RIV and the DRN package.
%
% The data for the 3D arrays are provided in separate file, one for each layer
% and contained in one or more zip files.
%
% Some special care had to be taken with respect to the elevations of the
% tops and bottoms of the layers. The top is contained in the AHN file, the
% "Algemene Hoogtekaart Nederland" (General Elevation MAP (or DEM) of the
% Netherlands). The elevations in it are given in m units above national
% datum (NAP), while all other layer-top and bottom elevations are given in cm.
%
% Anisotropy, necessary to model the ice-pushed formations in the center of the
% Netherlands is defined through an anisotropy factor and the angle of the
% direction of the main anisotropy axis with respect to the North
% (clockwise positive). Anioatropy may be neglected when modeling parts
% without it.
%
% The files do not contain basic data, but just the derived data
% values per model cell. There is one exception: The WEL list also
% contains the actual x,y coordinates of wells and their extraction.
% The extraction is some long-term average value (details to be looked up).
%
% The Excel workbook mflab/examples/NHI/NHI.xls lists all zip files and the
% ASCII files as well as their translation to the mat files
% generated from them. When the NHI model is renewed, the zip files have to be
% reloaded and this workbook may have to be updated as well. Currently, the
% entire model consists of about 120 ascii files.
%
% The NHI site and the NHI documentation on that site lack a descripition
% of all the individual files. This makes it hard to figure out what
% exactly is their contents, and even more so, how the files have been
% generated and how to interpret them correctly. It is the ambition to find
% out, so that their description can be added to the workbook HNI.xls for
% reference.
%
% Finally, note that the data provided by www.NHI.nu are only steady-state
% data. There are no dynamic stress period data available at this time.
% Exact 3D data on the subsurface salinity distribution are also lacking.
% I will do my best to use the 3D salinity data of the Deltares institute
% model for South-Holland, which covers about half the country.
% When the 3D data become available, it will be downloaded and will
% replace the South-Holland data.
%
% The current data set was downloaded around early April 2012.
%
% TO 120808

%% Variable types for layer processing
% Each varialble will be given a type that will be stored along with the
% variable in the matfile, to facilitate handling the variable when reloaded
% later on. The possible type so far are

varTypes = {
    '3Dlay';           % a 3D array size of model (Ny,Nx,NLay)
    '3Dcbd';           % 2 3D array size of confining beds (Ny,Nz,Ncbd)
    '3Dtime';          % a 2D-time array, e.g. recharge (Ny,Nx,Nt)
    'gridObj';         % a matlab grid object
    'wellObj';         % an array of mflab well objects
    'wellSeriesObj';   % an array of mflab well series objects
    'struct';          % a struct like ZETA
    'stress';          % a stres (GHB, RIV, ..) i.e. [it iz iy ix values]
    'zlist'            % one value per layer list such as wettability
    };

%% Setting up the directories and basename
% Basename used for this directory and the workbook with the file names
basename  = 'NHI';

%%
% Path to NHI dir, adapt for your personal computer
NHIdir  = '/Users/Theo/GRWMODELS/mflab/examples/NHI/';

asciiFilesDir  = 'NHIascii';    % name of subdirectory with NHI ascii files
datadir        = 'NHIdata';     % name of subdirectory with the NHI mat files
sheetNm        = 'files';       % name of worksheet in <<basename>>.xls with names of zipped and ascii files
columnHdr      = 'variable';    % name of column in <<sheetNm>> with MODFLOW variable name connected to ascii file

%%
% Inform user
fmt = '%56s = %s\n';
fprintf(fmt,'basename',basename);
fprintf(fmt,'subdirectory with NHI ascii files',NHIdir);
fprintf(fmt,'subdirectory with NHI mat files',datadir);
fprintf(fmt,'workshee with NHI ascii file names',sheetNm);
fprintf(fmt,'column naming the model variable or resulting mat file',columnHdr);

%% Change to NHIdir directory
eval(['cd ' NHIdir]);

%% Determine the number of NHI layers automatically
% The number of NHI layers equals the number of aquifers, taken from the files starting with kD
Nlay = length(dir(fullfile(asciiFilesDir,'KD*.ASC')));
fprintf(fmt,'The number of model layers',sprintf('%d',Nlay));

%% Determing the number of confining beds
% All NHI layers (except the last) have a confining bed below them, so we
% skip this by setting LAYCBD immediately as follows
LAYCBD = ones(1,Nlay);

%% Set minimum layer thickness
% Many layers have zero thickness over at least parts of their exent,
% whlile subtraction of the first bottom from the DEM even yields negative
% thickness over part of the area. Therefore, we guarantee a minimum
% thickness as specified by MINDZ here
MINDZ = 0.1;

fprintf(fmt,'LAYCBD','all but the last layers have a confining bed below them');
fprintf(fmt,'MINDZ: Minimum layer thickness used in NHI.\n',MINDZ);

%% Get meta data for cell-arrays
% Every ASCII file specifying some layer variable such as KD1.ASC contains
% meta data concerning the extent of the grid end its lower left coordinate
% in the national Dutch system ("Rijksdriehoeksmeting"). Note that conversion
% between the national x,y and the international WGS84 system used by Google
% can be done within mfLab using the WGS2RD and RD2WGS functions.
% Here we read out the meta data, using the ASCII file KD1.ASC for the
% purpose.

fid = fopen(fullfile(asciiFilesDir,'KD1.ASC'),'r');
    meta = getNHImeta(fid);
fclose(fid);

%% Compute grid-line coordinates for NHI based on meta data just obtained
% The cell center coordinates follow from the grid line coordinates and the
% cell size. We don't need them, because we will use mfLab's  grid object later
% on for all coordinate work.

fprintf('Computing the NHI grid line coordinates from NHI meta data.\n');

xGrNHI = meta.XLLCORNER+(0:+1:meta.NCOLS)*meta.CELLSIZE;
yGrNHI = meta.YLLCORNER+(meta.NROWS:-1:0)*meta.CELLSIZE;

%%
% Compute the indices of the model. Subsets may be used in functions
% applied further down to select data for a submodel.

Ix= 1:meta.NCOLS;  fprintf('Number of columns in NHI = %d\n',meta.NCOLS);
Iy= 1:meta.NROWS;  fprintf('Number of rows    in NHI = %d\n',meta.NROWS);

%% Get the actual data and save them into mat files in the data directory
% To generate the mat files for cell layer variables such as transmmissivity
% and elevations, we have to read out the ascii file of each layer and join
% them in a 3D MATLAB array. In the process we may rename the variable to
% a more international standard, like TRAN instead of kD, and we may
% convert the variable as required by MODFLOW, such as VCONT instead of C.
% Instead of separate 3D files of the tops and the bottoms of all layers,
% we use Z, which is a sequence of layer interface planes, starting at
% ground surface and ending at the bottom of the model, which combines both
% the aquifer-layers and the aquitard confining beds.
    
%% Get transmissivities
[TRAN,NHIvars.TRAN]  = makeArray(basename,sheetNm,asciiFilesDir,columnHdr,'TRAN',Ix,Iy);
TRANtype = '3Dlay';
save(fullfile(NHIdir,datadir,'TRAN'),'TRAN','TRANtype');

%% Vertical hydraulic resistance of confining beds
% In the Netherlands we always characterize the hydraulic vertical
% conductance combined with the thickness by their total hydraulic
% resistance, which equals c=DZ/kv (time]. In the USA and MOFLOW VCONT is used in
% stead, which is kv/DZ [1/time]. Therefore there is no fundamental difference
% between using one or the other. As the www.NHI.nu site publishes the
% c-values we will use them as well and take the reciprocal whenever we
% need them in MODFLOW later on.
%
[C NHIvars.C] = makeArray(basename,sheetNm,asciiFilesDir,columnHdr,'C',Ix,Iy');
Ctype = '3Dlay';
save(fullfile(NHIdir,datadir,'C'),'C','Ctype');

%% Starting heads
% Although not strictly necessary for steady-state models, unless used for
% fixed heads, we will still use the starting heads to initalize the the solver.
% This starts the convergence as closely as possible to the final results.
[STRTHD NHIvars.STRTHD]  = makeArray(basename,sheetNm,asciiFilesDir,columnHdr,'STRTHD',Ix,Iy);
STRTHDtype = '3Dlay';
save(fullfile(NHIdir,datadir,'STRTHD'),'STRTHD','STRTHDtype');

%% Cell-based recharge
% The NHI site gives one ASCII cell-based file for the net recharge. It
% must be some properly chosen long-term average value, which does not
% become clear from the site. The exact origin is still to be requested
% from the developers.
% The recharge is a 2D matrix. There is no need to generate a 3D array from
% it.
[fname NHIvars.RECH] = getNHIfileNm(basename,sheetNm,columnHdr,'RECH'); 
RECH     = getNHIASC(fullfile(asciiFilesDir,fname),Ix,Iy);
RECHtype = '3Dtime';

save(fullfile(NHIdir,datadir,'RECH'),'RECH','RECHtype');

%% Elevation of all layers and model grid object

%%
% Get layer top and bottom elevations of layers and aquitards
[Z NHIvars.Z] = makeArray(basename,sheetNm,asciiFilesDir,columnHdr,'Z',Ix,Iy);

%%
% Convert all Z(:,:,2:end) from cm to m, except Z1 = groundsurface which is
% already in m above datum
Z(:,:,2:end)=Z(:,:,2:end)/100;

%%
% Remove the NaNs of internal surface waters in AHN
% The AHN DEM contains nodata where surface water is present. For any model
% use we need proper elevation data here as well. We obtain these data by
% replacing the NaNs in the ASCII file by corresponding values of Z(:,:,2),
% the second layer, i.e. the bottom of the first model layer. Note that
% zero thickness that are the consequence of this procedure will be
% eleminated next.

z2=Z(:,:,2); % second plane

%%
% list of indices in the first layer correspond to second loose layer
Z(isnan(Z(:,:,1)))=z2(isnan(Z(:,:,1)));

%%
% Make sure that bottoms are at or below tops of layers working downward.
for i=2:size(Z,3)
    dz = max(diff(Z(:,:,[i i-1]),1,3),MINDZ);
    Z(:,:,i)=Z(:,:,i-1)-dz; 
end

%%
% Verify that each laye ris below its top. To do this, get DZ and and compute
% the minimum thickness of each layer.
DZ= abs(diff(Z,1,3));
for iz=1:size(DZ,3);
    fprintf('Thickness DZ of layer %2d varies between %12.3f and%12.3f\n',...
        iz,min(min(DZ(:,:,iz))),max(max(DZ(:,:,iz))));
end

%% 
% save Z
Ztype = '3Dlay';
save(fullfile(NHIdir,datadir,'Z'),'Z','Ztype');

%% Generate the mfLab grid object:
% The mfLab grid object contains all grid data and the methods to deal with
% them. It is a key object in mfLab. So create it:
gr= gridObj(xGrNHI,yGrNHI,Z,LAYCBD,MINDZ);
NHIvars.gr ={ 'gr' 'noFile' 'nozip' 'generated NHI grid' };

%%
% Save the grid
grtype = 'gridObj';
save(fullfile(NHIdir,datadir,'gr'),'gr','grtype');
    
%% Stresses
% Next, the stresses are read and stored. The stresses comprise the wells,
% and the data lists for the general head-boundary package (GHB), the drain
% package (DRN) the river package (RIV) and possibly other packages such as
% the drains with return flow (DRT), the changing head package (CHD) and
% the steam-routing package (STR). However, until now NHI only uses WEL, DRN,
% GHB and RIV packages.
% The stress packages are lists with one line per model cell. A list possibly
% contains hundreds of thousands of lines, each of which specifies the
% indices of the position of the cell in the grid and its data.
% Processing may, therefore, take same more time, although it is not generally
% overwhelming.
%
% When reading the stresses we have to consider the indices of the current grid,
% which are, in the case we select a submodel, a subset of all indices,
% requiring renumbering of the indices on every line of each of the stress lists.
% This renumbmering is done here to remain general. If the total model is
% selected, then renumbering has no effect.

%% Get the stresses for the RIV package (surface water)
% The NHI site gives several files for this package.Each refers to a different
% level of the represented surface water in the cells in a hierarchical
% surface-water ranking system. There are also files for different seasons,
% i.e. winter and summer and year-average files. We will only use the
% year-average files here at this time.
fprintf('Getting RIV for interaction with main, primary, secondary and tertiary surface water.\n');
[RIV0 NHIvars.RIV0] = getNHIBCN({'RIVHMY'  'RIVCM'  'RIVBMY' },basename,sheetNm,columnHdr,asciiFilesDir,Ix,Iy);
[RIV1 NHIvars.RIV1] = getNHIBCN({'RIVHP1Y' 'RIVCP1' 'RIVBP1Y'},basename,sheetNm,columnHdr,asciiFilesDir,Ix,Iy);
[RIV2 NHIvars.RIV2] = getNHIBCN({'RIVHS1Y' 'RIVCS1' 'RIVBS1Y'},basename,sheetNm,columnHdr,asciiFilesDir,Ix,Iy);
[RIV3 NHIvars.RIV3] = getNHIBCN({'RIVHT1Y' 'RIVCT1' 'RIVBT1Y'},basename,sheetNm,columnHdr,asciiFilesDir,Ix,Iy);

%%
% Put RIV into the layer that agrees with the corresponding depth of the cell bottom:
RIV0=setRIVdepth(RIV0,Z,gr.LAYCBD);
RIV1=setRIVdepth(RIV1,Z,gr.LAYCBD);
RIV2=setRIVdepth(RIV2,Z,gr.LAYCBD);
RIV3=setRIVdepth(RIV3,Z,gr.LAYCBD);

RIV0type = 'stress';
RIV1type = 'stress';
RIV2type = 'stress';
RIV3type = 'stress';
save(fullfile(NHIdir,datadir,'RIV0'),'RIV0','RIV0type');
save(fullfile(NHIdir,datadir,'RIV1'),'RIV1','RIV1type');
save(fullfile(NHIdir,datadir,'RIV2'),'RIV2','RIV2type');
save(fullfile(NHIdir,datadir,'RIV3'),'RIV3','RIV3type');

%% Get the stresses for the general head package (GHB)
% The structure of the files is as described before for the RIV
% package. The files hereafter are used to set the infiltration
% resistance. This is to differenciate infiltration resistance
% from drainage resistance.
% Therefore, RIV is used to extract the corresponding heads.
fprintf('Getting GHB to be used for infiltration resistances:\n');
[GHB0 NHIvars.GHB0] = getNHIBCN({'RIVHMY'  'GHBCM' },basename,sheetNm,columnHdr,asciiFilesDir,Ix,Iy);
[GHB1 NHIvars.GHB1] = getNHIBCN({'RIVHP1Y' 'GHBCP1'},basename,sheetNm,columnHdr,asciiFilesDir,Ix,Iy);
[GHB2 NHIvars.GHB2] = getNHIBCN({'RIVHS1Y' 'GHBCS1'},basename,sheetNm,columnHdr,asciiFilesDir,Ix,Iy);

GHB0type = 'stress';
GHB1type = 'stress';
GHB2type = 'stress';
save(fullfile(NHIdir,datadir,'GHB0'),'GHB0','GHB0type');
save(fullfile(NHIdir,datadir,'GHB1'),'GHB1','GHB1type');
save(fullfile(NHIdir,datadir,'GHB2'),'GHB2','GHB2type');

%% Get the stresses for the drain package (DRN)
% The NHI model uses the drain package to specify tile drainage only.
[DRN NHIvars.DRN]  = getNHIBCN({'DRNH' 'DRNC'},basename,sheetNm,columnHdr,asciiFilesDir,Ix,Iy);
DRNtype = 'stress';
save(fullfile(NHIdir,datadir,'DRN'),'DRN','DRNtype');

%% Get the wells (WEL)
% The NHI well file is special in that it contains two types of wells. The
% lines beyond about 3000 have as a comment
% the word 'SUR'. I, therefore, named these wells 'SURwells'.
% Plotting the surwells on a map reveals that they fully cover certain areas
% in the south-eastern part of the Netherlands. So presmumably these wells
% are irrigation wells, used to irrigate land from local groudwater.
% They are also shallow.
% The first about 3000 lines concern regulary wells, generally those
% of the drinking-water companies and industries. These wells are most often
% deeper, with screens in semi-confined aquifers. Each line of these wells
% contains further details, namely the actual well-specific data, among
% which are its x,y location and its actual extraction. More than one well may
% reside in a single cell, as MODFLOW will simply add them. This way the
% true well data are contained in the NHI well file.
%
% When reading the well file, a separate array of mflab well objects is
% generated with the actual
% well data. This array can be used with any model and submodel without the
% wel file list, as the latter is specific to the grid of the NHI groundwater
% model.
% Each of the lines in the regular wells also has a decrete value of the
% local transmmissivity. This too may be very useful for further modeling. A
% special transmissivy object array is also generated during the while
% reading and processing the well file.

%%
% First get well-file name
[fname names] = getNHIfileNm(basename,sheetNm,columnHdr,'WELL'); 

%% Get wells.
% Note that surwells(=irrigation wells?) are special, kD is at well sites.
% * well is the array of well objects
% * surwell is the array of surwell objects
% * kD is the array of transmissivit objects
% * WEL is the well list such as is processed by the MODFLOW's well
% package.
[well,surwell,kD,WEL] = getNHISCD(gr,fullfile(asciiFilesDir,fname),Ix,Iy);

NHIvars.WEL     = names;
NHIvars.well    = NHIvars.WEL;
NHIvars.surwell = NHIvars.WEL;
NHIvars.kD      = NHIvars.WEL;

welltype    = 'wellObj';
surwelltype = 'wellObj';
kDtype      = 'kDObj';
WELtype     = 'stress';
save(fullfile(NHIdir,datadir,'well'),'well','welltype');
save(fullfile(NHIdir,datadir,'surwell'),'surwell','surwelltype');
save(fullfile(NHIdir,datadir,'kD'),'kD','kDtype');
save(fullfile(NHIdir,datadir,'WEL'),'WEL','WELtype');

%% Done
fprintf('\n');
fprintf('All NHI data are saved as matfiles in the directory   <<%s>>\n',fullfile(NHIdir,datadir,''));
fprintf('File and variable names are in Excel workbook         <<%s>>\n',fullfile(NHIdir,[basename '.xls']));

%% Save the struct whit the names of the saved variables
NHIvarstype = 'vars';
save(fullfile(NHIdir,datadir,'NHIvars'),'NHIvars','NHIvarstype');

%% TODO:
% * Verification
% * Viewing
%
%% Show STRTHD in 3D
% For any part of the model coutour this part using indces and put this in
% a new grid (see how this in done in extactFromNHI.m in the example directories
% may be too big to show on smaller PCs

%h=showLayers(gr,STRTHD,well,'Countrywide model');
