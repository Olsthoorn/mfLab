%% Visualization of steady model with transport of a tracer and temperature

%% Close existing graphs and variables to prevent memory overflow
close all;
clear variables;

%% The options and steps will be outlined here
% We will use an animation to visualize the model results, followed by
% concentation and temperature breakthrough.
% Notice that there are many more visualization options.

%% Retrieve the basename of this model
load name        % retrieve basename stored in file name.mat
load(basename)   % get model arrays that were saved by mf_setup
load underneath  % info from mf_adapt, contains the species names

%% Generate an animation object
% The animation object created reads the simulated concentrations, the heads
% and budget (cell by cell flows). It needs the budget to draw the streamlines
% in the cross sections.
% Heads are only required if they are to be plotted.
animate = animateObj(basename,species,'head','budget');

%% Simulate the concentration of all species
animate.concXS(gr)

%% Generate observation wells and install them in the grid.
obsWells = observationObj(basename,'observations',gr,HK,[species,'head']);

%% Plot observations of given species. 
% The color, order of the graphs and their lineWidths are obtained from the
% worsheet observations" in workbook basename. This worksheet name was
% given in the call of observationObj above.

% Plot break-through graphs of species 1, all observation wells
obsWells.plot('fig',species{1},species{1},'ylabel',[species(1) ' [g/L]']);

% same for species(2), i.e. the temperature.
obsWells.plot('fig',species{2},species{2},'ylabel',[species(2) ' oC']);

% finished.

% TO 130614