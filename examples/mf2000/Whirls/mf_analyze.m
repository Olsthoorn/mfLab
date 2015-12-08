%% Analyzing output of the model
% TO 091011 091129
 
load('name.mat') % loads name.mat, contains the variable "basename"                 
load(basename);  % having retrieved baasename load the data in basename.mat
load('underneath'); % from mf_adapt

mf_checkdir;

H=readDat([basename,'','.hds']); % read the unformatted head file
hrange = ContourRange(H,50);

contourf(gr.xc,gr.yc,H(end).values(:,:,iM),hrange); % contour layer iLay

%% Showing particles

figure; hold on; view(3); xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('Particles (Whirls)');

gr.plotMesh('faceAlpha',0.15);

%% starting points
pGrp.plot('ro');

%% endpoints
pGrp = pGrp.getEndPoints('mp6.endp');
    %pGrp.dispEndPoints();
pGrp.endPointStatistics;
    %pGrp.plotEndPoints(1,'ko');

%% pathlines
pGrp = pGrp.getPathLines('mp6.plines');

pGrp.plotPath('b');

view(3);
