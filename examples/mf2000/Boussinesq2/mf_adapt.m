% Example Boussinesq (flow on an sloped base)
% TO 100414 100507 100917
%
% Boussinesq studied the flow of unconfined groundwater on a sloping aquifer base.
% This can be modeled with modflow by letting the cells step down.
%
% This file implements a simple one-dimensional unconfined model (along the
% x-axis, which steps down a given slope (the waviness of
% the slope can be adjusted with the parameter A. A=0 gives a straight slope.
% The mfLab userguide compares the results obtained by MODFLOW with the
% analytical solution given by Steward in WRR (2007).
%
% Note that the drying and rewetting mechanism of MODFLOW completely fails
% as it never converges. This is a long standing frustrating issue with
% respect to MODFLOW.
% Others (Doherty (2000) and Painter (2008) have developed better,
% more robust mechanisms, which are included in his MF2KASP version of
% MF2K (downloadable from www.pesthomepage.org).
%
% The current model defined below and in the accompanying workbook
% Boussinesq.xls, is transient, and uses Doherty's MF2KASP. It perfectly
% covnerges. Because of its usefulness I include MF2KASP with mfLab and
% would like to give full credit to Johns great work!
% 
% TO 100917

clear variables;

%% The model name. Every model has a "basename" and all associated files
basename='Boussinesq';

%% Parameters

k  = 10;                 % [md ] conductivity (uniform)
Sy = 0.1;               % [   ] specific yield
Ss = 0.0005;            % [m-1] elastic specific storage coefficient
z0 = 0;                 % [NAP] z at xm(1)
D  = 50;                % [ m ] dummy aquifer thickness (because unconfined)
slope = -1/10;          % [m/m] dB/dz, inclination with B the bottom of the aquifer
Q0 = 5;                 % [m2/d] discharge per m perpendicular to cross section

%% Specify grid line coordinates
dx  = 50;                                  % [ m ] grid cell width (uniform)
xGr = 0:dx:5000;                         % [m] cell coordinates
xm  = 0.5*(xGr(1:end-1)+xGr(2:end));    % [m] cell centers
yGr = [-0.5 0.5];                        % [m] cell coordinates (1 row model)

%% Aquifer bottom, top and center elevation

A  = 0;                 % [ m ] set A=0 for constant slope
zT = ones(size(xm))*D;  % [NAP] constant elevation of top of water table aquifer

% Parameterized elevation of wavy aquifer bottom of hill slope
zB = z0+slope*(xm-xm(1))+A*sin(2*pi*(xm-xm(1))/(xm(end)-xm(1))*2);  % Aquifer bottom

Nx = numel(xm);
Ny = 1;
Nz = 1;

Z = NaN(Ny,Nx,Nz+1);  % [NAP] elevation of model layer boundaries

Z(:,:,end)  = zB;
Z(:,:,1)    = zT+10*D;

gr = gridObj(xGr,yGr,Z);

%% Arrays
IBOUND = gr.const(1); IBOUND(:,[1 end],:)=-1; % Down-stream fixed head
HK     = gr.const(k);                     % [m/d] hor conductivity
VK     = gr.const(k);                     % [m/d] vertical conductivity
SS     = gr.const(Ss);                    % [1/m] elastic storage coefficient
SY     = gr.const(Sy);                    % [   ] specific yield

A=0;
STRTHD = gr.ZBlay+D + A* cos(pi*(gr.xm-mean(gr.xm))/diff(gr.xm([1 end])));                         % [ m ] initial head 1 m above aquifer bottom

well = wellObj(basename,'wells',gr,HK,{'PER','Q'});

%% Recharge
% see specified recharge in the workbook

[PERnams,PERvals] = getExcelData(basename,'PER','hor');
rch = mean(PERvals(:,strmatchi('RECH',PERnams)));

save underneath slope k Sy rch
