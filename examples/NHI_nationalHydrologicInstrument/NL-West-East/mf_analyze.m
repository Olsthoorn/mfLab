%% Analyzing output of the model (Part of NHI)
% This file analyzes (or visualizes) the results of the flow model that was
% invoked by typing mf_setup at the command line.
% 3D visualization is applied to allow inspecting complex 3D models like
% the one here.
% The model in his example is a 40 by 40 km submodel of the National Hyrological
% Instrument for the center of the Netherlands, around Lexmond.
% The data used have been downloaded from the official site www.NHI.nu/bibliotheek.
% The model is steady state, uses the average HI-provided recharge
% and also uses NHI-provided year-averaged boundary conditions (stresses).
% Summer and winter stresses are also available but have not been used
% given the lack of transient recharge data.
%
% TO 091011 120603 120816

%% Loading the basename, model data and extra info if available
load('name.mat') % get basename stored in file name.mat
load(basename);  % load the mat file with the model arrays stored in mf_adapt.m
load underneath  % load extra information if available

%% Print size and center coordinates of the model
fprintf('Size of model is %.0fkm in both directions\n',diff(gr.xGr([1 end])/1000));
fprintf('Center of model in Rijksdriehoek Coordinates (Dutch National system)\n');
fprintf(' is x=%.0f y=%.0f km\n',mean(gr.xGr)/1000, mean(gr.yGr)/1000);

[Lon,Lat]=rd2wgs(mean(gr.xGr),mean(gr.yGr)); %% In WGS coordinates
fprintf('In WGS (i.e. Google Earth):\nLon(Easting)=%g\nLat(Northing)=%g\n',Lon,Lat);

%% Load and show the computed heads of the model
H=maskH(readDat([basename,'.HDS']));

%%
% This shows the model layers in 3D. The figure starts with the heads which
% are H(1).values or H.values because the model is steady state.
% The 7 NHI layers are shown with a wire mesh around the entire model. This
% The wire mesh is drawn for 3 sides (west, south and top), but other sides
% can be added with the buttons on the figure. Each side can be switched on
% and off at will. The elevation of the top mesh corresponds to the
% elevation of ground surface and that of the bottom mesh with the bottom
% of the system (layer 7 in the current NHI).

%h=showLayers(gr,'H.values',well,basename);

%%
% By clicking on the rotation tool in the toolbar, one can turn the model
% in 3D space.
%
% To see the value of the colors, switch on the color bar in the toolbar.
%
% If wells are not well visible, change their color using the switch box
% near the right bottom of the screen.

%h=showLayers(gr,'H.values',well,basename,'contours');

%% Contour any values of the layers (head, transmissivities whatever)
% if exist('H','var'),    gr.layerContours('Heads [m]'                   ,H(end).values,basename,'m'); end

%%
% for widely varying positive values like transmissivties, take the log10
% if exist('TRAN','var'), gr.layerContours('log10(Transmissivity)'          ,log10(TRAN          ),basename,'log10(m2/d)'); end
% if exist('HK','var')  , gr.layerContours('log10(Transmissivity)'          ,log10(gr.DZlay.*HK  ),basename,'log10(m2/d)'); end
% if exist('VCONT','var'),gr.layerContours('log10(Vertical hydr resistance)',log10(1./VCONT      ),basename,'log10(d)'); end
% if exist('VKCB','var') ,gr.layerContours('log10(Vertical hydr resistance)',log10(gr.DZcbd./VKCB),basename,'log10(d)'); end

%% Show the specific vertical leakage downward in all layers

B   = readBud([basename '.BGT']); % read budget file with computed cell by cell flows
%%
% Compute specific vertical lekage as flowlowerface/cell area
% VL = B(end).term{strmatchi('FLOWLOWERFACE',B(end).label)}./gr.AREA3;
%%
% show it in contours
% gr.layerContours('Vert leakage [m/d]',VL(:,:,1:end-1),basename,'m/d');

%% Boundary conditions (stresses)
% Stresser are lists with one line per cell. To show then spatically they
% are put on a 3D cell grid and shown with one color in each cell. Taking
% column 5 of the lists implies the specified head or elevation of the stresses
% if exist('RIV','var'), gr.showStress('RIV',RIV,basename,1,5); end
% if exist('GHB','var'), gr.showStress('GHB',GHB,basename,1,5); end
% if exist('DRN','var'), gr.showStress('DRN',DRN,basename,1,5); end

%% Use zone budget to get budget overview
% - show the overall budget.
% - see help zonebudget for options to fine-tune
zonebudget(B);

%%

B=mf_Psi(B);

figure('position',get(0,'screensize')); grey = get(gcf,'color'); hold on;

gr.plotXSec(1,1:gr.Nlay,'hLines','k','smooth');

prange = ContourRange(B,25,[],'Psi');

plot(gr.xm,XS(H(1).values));

gr.streamlines(gca,B.Psi,prange,'color',grey);
