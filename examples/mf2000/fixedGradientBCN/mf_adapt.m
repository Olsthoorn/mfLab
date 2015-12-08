%% negConcductance
% experiment with negative conductances
% to generate fixed-head gradient boundary conditions.
% These fixed-head gradient boundary conditions are achieved using drains
% at the bottom of each boundary cell and computing the correct condutance
% for each of them.
%
% The prescribed gradient is obtained when the water table at the boundary
% is below the top of the aquifer and when the flow is outward.
%
% The method is developed for the Khettara model by Koert Strikker / Koen
% Geul for the Jorf area, near Erfoud, Morocco.
%
% TO 141104

clear variables; close all; % clean up

basename='fixedGrad'; % name of model and basename all its files

%% Generate the grid

% We use a cross section with 5 layers
xGr = 0:200:15000;
yGr = 0:1000:5000;
zGr = [200 -10 -20  -30 -40 -50];

gr = gridObj(xGr,yGr,zGr);

%% Aquitard and aquifer properties

kh   = [10; 20; 1; 5; 30]';  % layer conductances
%kh(:) = 30;
kv   = kh;                   % vertical layer conductivity
ss   = 0.0001;
sy   = 0.2;

%% Model arrays 3D
IBOUND = gr.const(1);  % Boundary array, >0 is which cells will be computed

STRTHD = gr.const(0);   % starting heads, all 0
STRTHD(:,1,:) = +30;

HK    = gr.const(kh);    % hor. conductivity
VK    = gr.const(kv);    % vertical conductivity
SS    = gr.const(ss);
SY    = gr.const(sy);

%% Set the boundary drains to achieve fixed gradient boundary conditions

zoneArray = gr.const(0);
zoneArray(:,1,:) = 1;  % left  drains
zoneArray(:,end,:)=2;  % right drains

K = diff(HK(:,:,end:-1:1),1,3);
K = cat(3,K(:,:,end:-1:1),HK(:,:,end));

gradientN = -0.6/1000;           % desired head gradient
gradientS = 1/1000;

hL = gr.ZBlay(zoneArray==1); % drain elevation left
hR = gr.ZBlay(zoneArray==2); % drain elevation right
C    = gr.DY(zoneArray==1)*gradientN.*K(zoneArray==1); % conductance drain

zoneVals = {1, hL, C };
%            2, hR, C};

% Generate drain array as required by modflow.
%DRN = bcnZone(basename,'DRN',zoneArray,zoneVals);

IBOUND(:,[1 end],:)=-1;

%north = lineObj(gr,basename,'boundaries','DRN','name','north');
south = lineObj(gr,basename,'boundaries','DRN','name','south');

%northG  = north.setGradient(gr,gradientN,HK);
southG = south.setGradient(gr,gradientS,HK);

clear north south

%% Get wells from worksheet 'PER'

well = wellObj(basename,'wells',gr,HK,'PER');

save underneath basename gradientN gradientS % zoneArray