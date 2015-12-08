% Example to show Toth's nested flow systems
% We start with the water table just hitting the lowest point in the area.
% The system gradually fills up, showing ever more complicatedd flow.
% TO 090806 091129 120424 130323

clear
clear variables; close all;

basename='Toth'; % model basename

kh   =  2.5;
kv   = 0.05;
peff = 0.35;
ss   = 1e-5;
sy   = 0.15;
Cdrn = 1000;
%% the grid

dx=100;
dz=10;

xGr= 0:dx:10000;
yGr= [-0.5 0.5];
zGr= 0:-dz:-250;

gr=gridObj(xGr,yGr,zGr);

%% zones and properties of subsurfac

IBOUND= gr.const(99);  % default, all 1 ==compute these cells
HK    = gr.const(kh);
VK    = gr.const(kv);
SS    = gr.const(ss);
SY    = gr.const(sy) * 0; % see below

% Because numerically the head is above the top of the digital model, the
% layers will never convert and Sy is inactive. One method to deal with
% this is to convert Ss of the top layer to Sy like Ss = Sy/Dz
SS(:,:,1) = sy./gr.DZ(:,:,1);

% RECH IS DEFINED IN WORKSHEET PER

%% generate an undulating landscape surface

% define ondulations: p [amplitude  waveLength[m] x_offset[m]]
p=[10  13000 500;  
   10    700  30;
   10    300 130];

h0 = 20;

% ondulating surface by adding an arbitrary number of waves as define in p
hdrn = h0 + sum(bsxfun( @times, p(:,1), sin( pi*bsxfun(@rdivide, bsxfun(@minus,gr.xm,p(:,3)) ,p(:,2))) ));

%% adapt model to landscape elevation

gr=gridObj(xGr,yGr,zGr); % regenerate the grid with adapted elevation

%% zone number for the drains --> IBOUND
idrn             = 10;
IBOUND(:,:,1)=idrn;

zoneVals = {idrn, hdrn, Cdrn};

%% show zones in IBOUND
%figure; showArray(gr,IBOUND,10,1,5);

%% Generate DRN for modflow
DRN = bcnZone(basename,'DRN',IBOUND,zoneVals);

%% Start at some value high enough to initially saturate the model
%  The drains at ground surface take care the flow will adapt to it.
STRTHD=gr.const(500);

save underneath hdrn