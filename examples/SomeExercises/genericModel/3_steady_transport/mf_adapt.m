%% CONCEPT Generic steady state transport model using MODFLOW2005 and MT3DMS

%% Explanation
%
% Notice that this model works, but needs more work to analyze it properly.
% The outcomes have not been analyzed yet to the extent that the model can
% be trusted. For now use the simple version in KonikowFig2 directory.
% TO 130618
%
% This model is the same as the steady generic transport model.
% modelling it has to be refined at least vertically. Furthermore, confining
% beds are no longer useful because the model needs to simulate also what happens
% inside aquitards. For this reason the layers in the worksheet LAY in the
% current workbook have been refined relative to the previous flow model. Calculation
% times will increase highly, especially if refining is also necessary
% areally. Finally we assume that the top of the aquifer is at 0 and that
% its top layers are not convertible (confined). We also change the drain
% elevation to zero, so that it is within the aquifer

close all;
clear variable;

basename = 'generic_steady_transport';
save('name','basename');

GREP = 'STRESS PERIOD';

xGr = 0:5000:75000;
yGr = 75000:-5000:0;

mf_lay2mdl(basename,xGr,yGr);

well = wellObj(basename,'wells',gr,HK,{'PER','Q'},{'PER','Tracer'});

IBOUND(:,1,[1 2]) = -1;

line = [ 7500 37500   0 ;
        47500 37500   0];

cDrn = 100; % [m/d]
DRN = gr.bcnLine(basename,'DRN',line,cDrn);


%% ======= TRANSPORT =========

%% Two species (constituents)
% In this example, we will use two species (dissolved consituents) in this example. The
% first one is conservative and the second is subject to sorption.
% Hence make sure NCOMP and MCOMP are both 2 in the MT3D worksheet.
% In the NAM sheet switch the packages BTN, ADV, DSP, SSM and GCG to ON.
% Because we will also simulate sorption, also switch RCT to ON.
% Because sorption will be linear, switch ISOTHM to 1 in the MT3D
% worksheet.
% For an accurate simulation choose the TVD solution method for the advection
% process. This is done by switching MXELM to -1 in the MT3D worksheet.
% Other parameters are given on a per layer basis in the worksheet LAY.
% Set dispersivity AL to 1.
% Also set bulk density RHOB = 1700.
% Set diffusion coefficients of both species DMCOEF = 1e-4 m2/d
% The reacton coefficients SP1 and SP2 are zero for species 1. The reaction
% coefficient SP1 for species 2 follows from linear sorption and a given 
% retardation of R5.
% Then
%    R = (\epsilon c + \rho_B K c)/(\epsilon c)
%      = 1+ \rho_B K \epsilon
% and so
%    K = \epsilon (R - 1)/ \rho_B
%      = .35*(5-1)/1700=8.235e-4
% This distribution coefficient has to be used for SP1, second species.

STCONC={};
STCONC{1} =  gr.const(1.0); % start conc of first constituent
STCONC{2} =  gr.const(1.0); % start conc of second consituent

ICBUND =  gr.const(1); % all cells will be computed

%save underneath zoneVals

%% Conclusion
% Particle tracking can be added to a model in just a few steps.
%% TO 130614