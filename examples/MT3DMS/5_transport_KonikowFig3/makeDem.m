%% Generate suitable top and base of this synthetic model
% The model was given by Konikow (2011), Groundwater Journal.

%% Explanation
%
% This file generates the contours of the lake, obstacles, river and model
% outline given in fig3 of Konikow (2011). These contours have been
% digitized on screen from the image 'KonikowFig3.png' and stored in mfiles
% lake.mat, Mdl.mat, Isl1.mat, Isl2.mat and Riv.mat.
% These files should not be deleted.

% A suitable top and bottom of the aquifer is computed as follows. First a
% confined model was made with the lake and the river as boundary
% conditions, and IBOUND set to 1, so that head contours are obtained
% everywhere in the image.

% The obtaind heads are used as top of the aquifer and its bottom is set to
% 15 m below. This Z array is stored in Z.mat and used to generate the
% gridObj in mf_adapt.

% The well and observation well depths have been obtained from this Z
% array, making sure the well screens are inside the aquifer.

close all;
clear variables;

basename = 'transportKonikowFig3';
save('name','basename');

GREP = 'STRESS PERIOD';

xGr = 0:100:4500;
yGr = 5400:-100:0;

mf_lay2mdl(basename,xGr,yGr); % generates model input arrays

% we use the initial grid with Z is [0 -15]. This is ok because we run the
% model in confined mode.

STRTHD(:) = 100;  % immaterial, because boundaries are set by lake and riv

load lines; % constains xLake yLake xMdl yMdl xIsl1 yIsl1 xIsl2 yIsl2 xRiv yRiv

% We only need the lake and the river to set head boundary conditions
iLake = 2; lakeStage = 75;                  % lake stage = +75 m
iRiv  = 3; rivStage = [2 26]; cRiv = 1000;  % rive stage from 2 to 26 m

% Draw the river line through the grid
P     = gr.lineObjects([xRiv,yRiv,-7.5*ones(size(yRiv))]);

L     = cumsum([P.L]');  % total river length
zRiv  = interp1(L([1 end]),rivStage,L);  % stage along river
river = [[P.xm]' [P.ym]' zRiv];          % 3D river coordinates

[iidx,I]=unique([P.idx]);  % remove double indices
river = river(I,:);        % remove double indices

%% Set IBOUND
IBOUND(inpolygon(gr.XM,gr.YM,xLake,yLake)) = iLake;
IBOUND([P.idx])                            = iRiv;

rivDepth = 5; 
chdVals = {iLake  lakeStage lakeStage};                     % head values for CHD
rivVals = {iRiv   river(:,end) cRiv river(:,end)-rivDepth}; % river input for RIV

CHD = bcnZone(basename,'CHD',IBOUND,chdVals);  % input for CHD
RIV = bcnZone(basename,'RIV',IBOUND,rivVals);  % input for RIV

return
%% After the model has run successfully

H = readDat([basename '.HDS']); % read the computed heads

D = 15; % Depth of aquifer
Z = H(end).values            % compute top elevation of new aquifer
Z(:,:,2) = Z(:,:,1)-D;      % compute bottom elevation of new aquifer

save dem.mat Z  % save Z as dem.mat