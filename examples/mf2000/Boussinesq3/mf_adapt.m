% Example Boussinesq (flow on an sloped base)
% TO 100414 100507 100917
%
% Boussinesq studied the flow of unconfined groundwater on a sloping aquifer base.
% This can be modeled with modflow by letting the cells step down.
%
% This file implements a two-dimensional unconfined model that delinces along the x-axis
% 500 m over 5000 m
% At the center there s a 100x100 m area with 2 m/d recharge.
% The area is drained. That is, drains exist along the top and downside boundaries
% This problem was proposed by Kick Hemker if the Boussinesq discussion in the
% Modflow user group. 
% TO 100927

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

iDrn = 5;

kStar = k*cos(slope)^2;

%% Specify grid line coordinates
xGr = (0:+50:5000);                         % [m] cell coordinates
yGr = (5000:-50:0)';                        % [m] cell coordinates 
[xGr,yGr,xm,ym,dx,dy,Nx,Ny] = modelsize(xGr,yGr);

Nz = 1;

%% Aquifer bottom, top and center elevation
zT = ones(Ny,Nx)*D;  % [NAP] constant elevation of top of water table aquifer

% Parameterized elevation of wavy aquifer bottom of hill slope
zB = z0+slope*(xm-xm(1));
zB = bsxfun(@times,ones(size(ym)),zB);


Z = NaN(Ny,Nx,Nz+1);  % [NAP] elevation of model layer boundaries

Z(:,:,end)  = zB;
Z(:,:,1)    = zT+10*D;

gr = gridObj(xGr,yGr,Z);

%% Arrays
IBOUND = gr.const(1); IBOUND(:,[1 end],:)=-1; % Down-stream fixed head
HK     = gr.const(kStar);                     % [m/d] hor conductivity
VK     = gr.const(kStar);                     % [m/d] vertical conductivity
SS     = gr.const(Ss);                    % [1/m] elastic storage coefficient
SY     = gr.const(Sy);                    % [   ] specific yield

STRTHD = NaN(Ny,Nx,Nz);
STRTHD(:,:,  1) = gr.ZBlay+D;
STRTHD(:,:,end) = STRTHD(:,:,1);
well = wellObj(basename,'wells',gr,HK,{'PER','Q'});

IBOUND(:,[1 end],end) = iDrn;

Cdrn = 1000;
DRN = bcnZone(basename,'DRN',IBOUND,{iDrn gr.ZBlay(IBOUND==iDrn),Cdrn});

%% Recharge
% see specified recharge in the workbook

rch = 2;
NPER = getPeriods(basename);
RECH = zeros(Ny,Nx);
RECH(gr.ym>2450 & gr.ym<2550, gr.xm>2450 & gr.xm<2550) = rch; 
RECH = bsxfun(@times,RECH,ones(1,1,NPER));

figure; spy(RECH);

save underneath slope k Sy rch
