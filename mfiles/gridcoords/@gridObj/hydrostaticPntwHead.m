function [STRTHD,p] = hydrostaticPntwHead(o,STRTHD,STCONC,drhodc,rho0)
%% [STRTHD,p] = gr.hydrostaticPntwHead(STRTHD,STCONC,drhodc,rho0);
% Only STRTHD(:,:,1) is used from the input
% Compute hydrostatic pointwater head whre you desire no vertical flows
% due to density variations. (i.e. if you want environmental heads).
% This is equivalent ot using CHD with CHDDENSOPT==2.
% It computes the pressure based on the head at the top of the STRTHD array
% and the given density only, taking into account the thickness of the
% layers and the head of the top layer above the first cell. Contrary to
% Seawat, IBOUND is ignored, so inactive cells are taken into account in
% this computation of the hydrostatic head. In those cases, the STCONC must
% be specified also for inactive cells.
% Only the top head of STRTHD is used and unaltered.
%
% TO 120628 120704

g = 10; % infact a dummy, don't need it here

if nargin<5, rho0=1000; end

rho = rho0 + drhodc*STCONC;

rho_g_dz = NaN(o.size);

% integral of rho g z from head to center of first layer
rho_g_dz(:,:,1    ) =  g * rho(:,:,1) .* (STRTHD(:,:,1)-o.ZMlay(:,:,1));

% additional from each layer center to the next layer center
rho_g_dz(:,:,2:end) =  g *(rho(:,:,1:end-1).*o.DZlay(:,:,1:end-1)+rho(:,:,2:end).*o.DZlay(:,:,2:end))/2;

% hydrostatic head
p = cumsum(rho_g_dz,3);

% STRTHD computed from hydrostatic pressure as point-water head so that
% Seawat can use it immediately. It overrides the initial STRTHD for all
% layers except the top.
STRTHD = o.ZMlay + p./(rho*g);
