%% NHI submodel Haarlemmermeer (Haarlem lake polder (Schiphol area))
% This file models an area around the location near Lijnden, Haarlem lake polder.
% The file generates a finite difference model from the date of the
% national groundwater model NHI as published on www.NHI.nu.

clear; close all;

grey = [0.8 0.8 0.8];
mustShowXSections=1;

%% Center the model grid using a placeholder set in Google Earth
% The point has been obtained by selecting a ping in Google Earth and
% saving it to a kml file. This kml file is then copied to the current
% directory
d = dir('*.kml'); % get file properties

[Path Location Ext]= fileparts(d(1).name); % split its name

%% Use the name of the kml file basename as the name of the model
basename = Location;

%% Get the wgs84 coordinates of the location
[E, N] = kmlpath([Location Ext]);  % wgs84 (Google Earth coordinates)
[xw,yw] = wgs2rd(E,N);             % Dutch national system coordinates (Rijksdriehoek)

%% Round center coordinates to center of cell in NHI grid

NHIcellsize = 250; % [m]

xC = NHIcellsize * floor( xw/NHIcellsize ) + NHIcellsize/2;
yC = NHIcellsize * floor( yw/NHIcellsize ) + NHIcellsize/2;

%% Choose size of model grid in terms of number of cells
ncells = 40; % means that the model will be 2*ND*NHIcelsize

xlim = xC + [-ncells ncells] * NHIcellsize;
ylim = yC + [-ncells ncells] * NHIcellsize;

%% Extract the data arrays for the submodel from the NHI arrays

Model = modelObj(['..' filesep 'NHIdata'],xlim,ylim);
Model = Model.descr(sprintf('Generated: Model %s as a submodel from NHI, %s',basename,datestr(now)));
Model = Model.BCF2LPF();
Model = Model.descr('Converting: BCF --> LPF');

Model.showBox(Model(1).description{end});

%% unpack the model object array to the workspace
unpack;

%% Show cross sections with only the aquitards colored

if mustShowXSections
    ND   = 20;  %distance between cross section plots in cells
%% Plotting section along x-axis
    IyXS= 1:ND:gr.Ny;
    gr.plotXSec(IyXS,[],'hlines',grey,'title',basename);
%% Plotting section along y axis
    IxYS= gr.Nx:-ND:1;
    gr.plotYSec(IxYS,[],'hlines',grey,'title',basename);   
end

%% Add information that is missing in NHI, because NHI data are only provided for cells
% These values were not be set during creation of the objects, because of
% the large time to compute for the entire NHI model (several minutes)
for i=length(well):-1:1,
    well(i).idx  = cellIndex(well(i).ix,well(i).iy,well(i).iLay,gr.size);
    kD  (i).idx  = well(i).idx;
    well(i).ztop = gr.ZTlay(well(i).iy,well(i).ix,1);
    well(i).z    = [gr.ZBlay(well(i).idx),gr.ZTlay(well(i).idx)];
    kD(  i).z    = well(i).z;
    well(i).DZ   = abs(diff(well(i).z));
    kd(  i).DZ   = well(i).DZ;
end

%% Select subset of all stresses for this model
RIV = RIV0;   clear RIV0 RIV1 RIV2 RIV3
GHB = GHB0;   clear GHB0 GHB1 GHB2 

%% Show STRTHD in 3D
% h=showLayers(gr,'STRTHD',well,Location);

%% Set the layer to 2 for all RIV, GHB and DRN that have a layer thickness < 1 m
%  This is necesary make rewetting work for thin layer above the water table.    
RIV = gr.mvStressDown('RIV',RIV);
GHB = gr.mvStressDown('GHB',GHB);
DRN = gr.mvStressDown('DRN',DRN);
    
%% Get IBOUND
IBOUND=gr.const(1);
IBOUND([1 end],:,:) = -1;
IBOUND(:,[1 end],:) = -1;
IBOUND(isnan(STRTHD))=0;

save underneath Location kD
