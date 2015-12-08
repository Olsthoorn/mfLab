%% Visualization of steady transport model Konikow (2011) Fig3, Groundwater Journal

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
animate = animateObj(basename,species,'head','budget');

%% Simulate the concentration of all species
animate.concXY(gr,{'Tracer','head'},1,well,'backgr','bckgrd','contourClr','y')

%H = readDat([basename '.HDS']);
%contour(gr.xm,gr.ym,H(end).values,0:5:80,'k');

%% Generate observation wells and install them in the grid.
obsWells = observationObj(basename,'observations',gr,HK,[species,'head']);

%% Plot break-through graphs of species 1, all observation wells
obsWells.plot('fig',species{1},species{1},'ylabel',[species(1) ' [g/L]']);

%% Finished.

% TO 130614