%% Horizontal well
% Example horizontal well, to check the capacity of a horizontal
% well drilled and installed in Nieuwegein, The Netherlands in March 2010
% using horizontal drilling techniques. These wells are called HHDW,
% Horizontal Direction Drilled Wells. The well was installed by consortium
% with the help of a research grant by the Senter Novem the Dutch ministry
% of Economic Affairs
%
% The 30 m long finished yielded 19 m3/h at 3.5 m drawdown.
% Further development should enhance its capacity, hopefully to its
% theoretical value, but was never realized due to cost.
%
% This mfLab example implements a model to a priori estimate the ideal
% capacity achievable when all resistance due to drilling mud would be removed.
%
% The well 5 m below the local confining bed with no leakage.
% Fixed head boundaries were set at 1000 m distance.
%
% We may compare the result with the analytical solution of
% Bruggeman, 1999, p364m solution 522.04.
%
% HDDW.m implements Bruggeman's analytical solution for the
% steady state situation.
% mf_adapt implements the same problem in MODFLOW through mfLab.
%
% Both solution use the same 3D grid to evaluate the drawdown, which
% facilitates comparison.
%
% The mfLab model is a rectangular grid of 1000 by 1000 m and 50 m deep.
% The grid is refined where required to more accurately model de detailed
% flow in the vicinity of the HDDW.
%
% The model grid is a quarter because of the symmetry.
%
% As always, you first run mf_setup and when the model finishes normally,
% run mf_analyze. mf_setup runs the local mf_adapt in which the model
% arrays are assembled.
% You can always run mf_adapt separately to verify that it works correctly.
% When ok, run mf_setup.
%
% TO 100610 130322

basename='HDDW';

k    = 15;         % [m/d]  hydraulic conductivity of the aquifer
L    = 17.5;       % [m]    half the length of the HDDW
peff = 0.35;       % [-]    effective porosiyy

% Model grid

yGr= sinespace(-1000,1000,70, -pi/2,pi/2);
xGr= sinespace(    0,1000,70,     0,pi/2);
zGr= [0 -0.1 -4.9 sinespace(-5, -50, 25, -pi/10,pi/2)];

gr = gridObj(xGr,yGr,zGr);

IBOUND = gr.const(1); IBOUND(:,:,1)=-1; IBOUND(:,end,:)=-1;
HK     = gr.const(k);
VK     = gr.const(k);  VK(gr.ZM>-5)=k/100;
PEFF   = gr.const(peff);
STRTHD = gr.const(0);

well = wellObj(basename,'wells',gr,HK,'PER');

save underneath peff  % Remember this for the file mf_analyze.m
