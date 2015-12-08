% show results of CFP example
%
% Read out the head and budget file and show the results
% This file is still a bit primitive and geared to the steady-state
% solution. It can however, be easily extended to other more extensive
% situations. the file is a part of the mflab construction and uses the
% modflow as the conduit low package
%
% TO 101231

clear variables; close all;

load ('name.mat');
load(basename);
load underneath
 

%% Get heads of first layers

H=readDat([basename '.HDS']);

%% Compute budgets
B=readBud([basename '.BGT']);

%% heads plot

% iy=DRN(1,3);  % where is the drain in the grid?

figure; hold on;

[~,hdl]=contourf(gr.xc,gr.yc,H(1).values(:,:,1),50,'edgecolor','none');

xlabel('x [m]'); ylabel('y [m]'); title('head contours qanat system');
colorbar;

figure; hold on

plot(gr.xm, gr.Z(iy,:,end),'k');
plot(gr.xm, gr.Z(iy,:,1)  ,'k');
patch([gr.xm gr.xm(end:-1:1)], [gr.Z(iy,:,end) gr.Z(iy,end:-1:1,1)],[1 0.9 0]);

plot(gr.xm,H(1).values(iy,:,1));
plot([P.xm],[P.zm],'r');

zonebudget(B);
