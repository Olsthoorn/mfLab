%% Generic steady state MODFLOW model with transport using MT3DMS

%% Explanation
%
% This model mimics the example given by Konikow (2011), fig 3.
% This model simulates a pollution in a groundwater field with a fair
% hydraulic gradient. Konikow uses it to show differences in outcomes
% caused by different solver methods.
% The model is single layer, 4.5 km wide and 5.4 km long with uniform cells
% of 100x100 m. Flow is from a high-level lake in the north along some
% obstacles to a river in the south. A pollution is simulated by means of
% two injection wells.

%% Steps described already in the generic steady model
close all;
clear variable;

basename = 'transportKonikowFig3';
save('name','basename');

GREP = 'STRESS PERIOD';

xGr = 0:100:4500;
yGr = 5400:-100:0;

mf_lay2mdl(basename,xGr,yGr);

load dem;  % contains Z (top and bottom of system
load lines % contains xLake yLake xMdl yMdl xIsl1 yIsl1 xIsl2 yIsl2 xRiv yRiv

gr = gridObj(xGr,yGr,Z);

iLake = 2; lakeStage = 75;
iRiv  = 3; rivStage = [2 26]; cRiv = 1000;

P     = gr.lineObjects([xRiv,yRiv],1);  %,-7.5*ones(size(yRiv))]);
river = [[P.xm]' [P.ym]' [P.zm]'];
[~,I]=unique([P.idx]);
river = river(I,:);

STRTHD(:) = 100;

IBOUND = gr.const(0);
IBOUND(inpolygon(gr.XM,gr.YM,xMdl,yMdl))   = 1;
IBOUND(inpolygon(gr.XM,gr.YM,xLake,yLake)) = iLake;
IBOUND(inpolygon(gr.XM,gr.YM,xIsl1,yIsl1)) = 0;
IBOUND(inpolygon(gr.XM,gr.YM,xIsl2,yIsl2)) = 0;
IBOUND([P.idx])                            = iRiv;

chdVals = {iLake lakeStage lakeStage};
rivVals = {iRiv river(:,end) cRiv river(:,end)-5};

CHD = bcnZone(basename,'CHD',IBOUND,chdVals);
RIV = bcnZone(basename,'RIV',IBOUND,rivVals);

species = {'Tracer'};

well = wellObj(basename,'wells',gr,HK,{'PER','Q'},species);

%% Transport modeling of two species (constituents)
%
% In this example, we will use two species (dissolved consituents). The
% first one is conservative and the second is subject to sorption. We
% choose its parameters such that it behaves like heat transfer. Therefore,
% we have a tracer and temprature. The retardation of temperature is about
% 2.
% Make sure NCOMP and MCOMP are both 2 in the MT3D worksheet.
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

STCONC{1} =  gr.const(0); % start conc of first constituent
STCONC{2} =  gr.const(0); % start conc of second consituent

ICBUND =  gr.const(1); % all cells will be computed

save underneath species  % store species names, to use it in mf_analyze

% TO 130614