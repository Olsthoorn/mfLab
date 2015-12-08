%% Discription 
% example of the Conduit Flow Process (CFP)
% By Amir Haidari AH (student at TU-Delft)
% control and corrected by: prof. Theo Olthoorn (profssor at TU-Delft)
% supported by: Waternet Amsterdam (The Netherlands)
% date: 110513 
%
% This example is from the CFP manual (page 33).
%

clc; clear variables; close all;
basename='exCFP';

%% Basic grids and basic parameters of layers 
LX   = 800;  % [ m ] lengh of aquifer 
LY   = 800;  % [ m ] width of aquifer 
ZT   = 7.0;  % [mNAP] elevation of ground surface at outlet of drain
ZB   =-0.2;  % [mNAP] elevation of bedrock at outlet of drain
DX   = 200;  % [ m ] cell size in x-direction
DY   = 200;  % [ m ] cell size in y-direction
Nz   =   3;  % [ - ] number of model layers
k    =  10;  % [m/d] conductivity (layers)
peff = 0.35; % [ - ] effective porosity (layers)
sy   = 0.25; % [ - ] specific yield (layers)
ss   = 1e-5; % [1/m] specific elastic storage (layers)

%% Cell size geometry apprach
%  See excel sheet 'CFPN'

xGr=0:DX:LX;  % [-] determines x boundary of each cell
yGr=0:DY:LY;  % [-] determines y boundary of each cell     

zGr=ZT:-(ZT-ZB)/Nz:ZB;

gr = gridObj(xGr,yGr,zGr);

%% IBOUND Variable
IBOUND = gr.const(1); IBOUND(end,:,:)=-1;
STRTHD = gr.const(0); STRTHD(IBOUND==-1)=7;

HK   = gr.const(k);
VK   = gr.const(k);
PEFF = gr.const(peff);
SY   = gr.const(sy);
SS   = gr.const(ss);

%% Show the model

figure; 
ax = axes('nextplot','add','xLim',gr.xGr([1 end]),'ylim',gr.yGr([end 1]),'zlim',gr.zGr([end 1]));
title(sprintf('My first attempt %s\n',datestr(now)));

gr.plotMesh('faceAlpha',0.15);
view(3)

         %% GET CFP data for plotting
[CFPNnams,CFPNvals] = getExcelData(basename,'CFPNODES','H');
[CFPPnams,CFPPvals] = getExcelData(basename,'CFPPIPES','H');

for ip=1:size(CFPPvals,1)
    I=CFPPvals(ip,[2 3]);
    plot3(gr.xm(CFPNvals(I,2)),gr.ym(CFPNvals(I,3)),squeeze(gr.zm(CFPNvals(I,4))),'r');
end

for in=1:size(CFPNvals,1)
    plot3(gr.xm(CFPNvals(in,2)),gr.ym(CFPNvals(in,3)),gr.zm(CFPNvals(in,4)),'bo');
end

%% Recharge only in the cells where water entrers pipe nodes

[pernams,pervals,NPER]=getPeriods(basename);
rech = gr.const(pervals(:,strmatchi('RECH',pernams)));

RECH=zeros(size(IBOUND(:,:,1)));
RECH(gr.Ny*(CFPNvals(:,2)-1)+CFPNvals(:,3))=1.0*(CFPNvals(:,8)>0);
RECH=repmat(RECH,[1,1,NPER]);

for iper=1:NPER, RECH(:,:,iper)=RECH(:,:,iper)*rech(iper); end


save underneath CFPNnams % anything you need in mf_analyze