%% NHI submodel NL-West-East
% Cross section model running W-E through the center of the Netherlands
% between the North Sea and the border with Germany.
%
% The model is a submodel extracted from the from the NHI, the Netherlands
% Hydrologic Instrument.
%
% The model is set up and run in the open-source mfLab environment, see
% http://code.google.com/p/mfLab

%% About the model
% The West-East cross section through the Netherlands is a famous one, as it has
% been published already in hydrogeological geological cross sections
% included with the geological maps of the Netherlands in the 1970s. Thanks
% to the new NHI, it is now possible to model it in its entirety, bringing
% it to life and allowing to draw stream lines in it to demonstrate the
% actual flow systems, which, until now had to be skeched based on
% hydrological insight alone.
%
% The cross section of the Netherlands has also been published by Dufoour
% (1998, 2000), stressing the interaction with fresh water floating on salt
% water throughout the Netherlands and how this interface waves due to the
% groundwater flows and heads as it originates from the landscape. This too
% can now, finally, be modelled using SEAWAT and or SWI (salt water
% intrusion package freely available on the internet).
%
% The grid used for the NHI is unsuitalbe for transport modeling with SEAWAT,
% for computation of displacement of salt and brackish water with dispersion,
% as it has only a limited number of layers, a number of which with zero
% thickness over at least part of their extent. We therefore adapt the
% standard NHI grid layer elevations to make it suitable for transport
% modeling, whithout changing transmissities and hydraulic resistances as
% we go.
%
% The first adaptation necessary for transport modelling is to change
% confining beds into regular model layers. MODFLOW discerns layers from
% so-called confining beds.  Models are layers in which the heads and the
% flows across cell faces are computed. Confining beds are layers without
% this information. In MODFLOW every layer may have a confining bed
% associated with it. The confining bed is connected to the bottom of the
% Layer. In fact, it consists only of a vertical resistance, between two
% vertically stacked layers, which is enough to compute heads in all layers
% and the vertical flows through the confining beds between each adjacent pair
% of model layers. This setup is used in the NHI, which has 7 model layers
% and 6 confining beds in between. This setup is unsuitable for transport
% modeling. The first adaptation is to change the confining beds into
% regular model layers, so that transport can also be computed within these
% confining layers. This changes the model in one with 13 model layers.
%
% In subsequent steps we will adapt the thickness of the model layers, to
% make sure that each layer has a convenient mininum thickness and that
% the layers are sufficiently refined vertically to make them suitable for
% transport modeling. The different steps will be shown graphically.

clear; close all;

%% Center the model grid the location of the pumping test.
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

%% Define grid

load(['..' filesep 'NHIdata' filesep 'gr.mat']);

CELLSIZE = gr.dx(1);

xShore = 84750;  % for this cross section to prevent NaNs
xLim = xC + [-100000 90000]; xLim(1) = xShore; 
yLim = yC + CELLSIZE *[-1 1];

clear gr

%% Extract the data arrays for the submodel from the NHI arrays

Model = modelObj(['..' filesep 'NHIdata'],xLim,yLim);
Model = Model.descr(sprintf('Generated: Model %s as a submodel from NHI, %s',basename,datestr(now)));
ic = strmatchi('C',{Model.name},'exact');
Model(ic).name = 'VCONT';
Model(ic).var  = 1./Model(ic).var;

%Model = Model.BCF2LPF();
%Model = Model.descr('Converting: BCF --> LPF');
%Model.showBox(Model(1).description{end});

%% unpack the model object array to the workspace
unpack; % also contains the grid

IBOUND=gr.const(1);

HY = TRAN./gr.DZlay;

%% Show cross sections with only the aquitards colored
grey = [0.8 0.8 0.8];
fprintf('Plotting cross sections along x-axis ...\n');
IyXS= hit(gr.yGr,yC);

gr.plotXSec(IyXS,'hlines',grey,'title',basename,'smooth');

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

if 1   % set 1 to 0 to demonstrate that this action is necessary
    %% Set the layer to 2 for all RIV, GHB and DRN that have a layer thickness < 1 m
    %  This is necesary make rewetting work for thin layer above the water table.    
    RIV = gr.mvStressDown('RIV',RIV);
    GHB = gr.mvStressDown('GHB',GHB);
    DRN = gr.mvStressDown('DRN',DRN);
end

%% SWI specific matrices
SSZ     = gr.const(0.35);        % effective porosity called SSZ in SWI
ISOURCE = gr.const(  -2);   % source zone incase of sinks 

%% ZETA planes

zta = [
     83746     -18.19
     85896     -76.07
     88226     -88.48
     91272     -95.37
     95932     -77.45
    100950     -40.24
    108656     -16.81
    126756     -25.08
    138226     -69.18
    146290    -113.97
    154534    -115.35
    164928    -105.70
    166900     -27.14
    172276     -19.56
    179982     -88.48
    188405    -144.98
    195036    -140.85
    197366     -87.10
    202921     -31.28
    211344     -10.60
    219946     -17.50
    227832     -19.56
    235896      -0.96
    241631      12.83
    248441      17.65
    255789      27.30];

%%
ZETA = gr.zeta(interp1(zta(:,1),zta(:,2),gr.xm));

[Nnum,Ntxt]=xlsread(basename,'NAM','','basic');
if Nnum(strmatchi('SWI',Ntxt(:,1),'exact'),3)
    SWIon = true;
else
    SWIon = false;
end


save underneath Location kD SWIon
