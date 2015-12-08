%% Qanat under construction
% TO 100827

%% Read out the head and budget file and show the results
% This file is still a bit primitive and geared to the steady-state
% solution. It can however, be easily extended to other more extensive
% situations.

clear variables; close all;

load name
load(basename)
load underneath
 
%% Get heads of first layers 

H = readDat([basename '.hds']);
H = maskHC(H,[-Inf 1000],[NaN NaN]); % get heads for these dates only and only layer 
B = readBud([basename '.bgt']); % selected cross sections for water balance

hrange = ContourRange(H,50);

%% heads plot
figure; hold on;
xlabel('x [m]'); ylabel('y [m]');
title('Head contours qanat system');

contourf(gr.xc,gr.yc,H(end).values(:,:,end),hrange);
plot(gr.XM(DRNidx),gr.YM(DRNidx),'b','lineWidth',3);

%%

%% Draw cross section showing model layers

gr.plotXSec(hit(gr.yGr,0),'fig','all','smooth','faceAlpha',0.25,'edgecolor','none');

xlabel('x [m]'); ylabel('elev [m]');
title('head through the qanat (at bottom of aquifer)');

plot(gr.XM(DRNidx),gr.ZM(DRNidx),'k--','lineWidth',1);

% Head along qanat
plot(gr.xm,H(end).values(hit(gr.yGr,0),:,end),'b','linewidth',2);
plot(gr.xm,H(end).values(      1      ,:,end),'b','linewidth',2);

% Draw the qanat
plot(xQanat,zQanat,'r','linewidth',3);

%% Zone budget is probably enough to get an overall budget over all items                           

zonebudget(B);
