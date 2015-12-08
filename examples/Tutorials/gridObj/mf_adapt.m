%% Construction of Tafilalt model Coert Stikker Juni 2013

fprintf('\n\n\n\n***** %s -- started *****\n\n',mfilename); tic;

%% Name of this project
basename= 'testCase';
save ('name','basename');

testData = 'testCaseData'; 

minDZ = 0.1;          % min allowed layer thickness

hk    = [0.4,2,100];  % [m/d] hor conductivity of the 3 geological ESRI layers
sy    = 0.25;         % [ - ] specific yield for transient simulations
ss    = 1e-5;         % [1/m] elastic storativity

xGr = 0:10:1000;
yGr = 1000:-10:0;
zGr = [0 -25 -50 -75];

gr = gridObj(xGr,yGr,zGr,'minDZ',minDZ);

DEM = gr.ZTlay(:,:,1) + 10*cos(3*gr.Xm/(1000*pi)).*sin(3*gr.Ym/(1000*pi)); %DEM


figure('name',basename,'position',screenPos(0.75));
ax = axes('nextPlot','add');

[~,h] = contourf(ax,gr.xc,gr.yc,DEM,30,'edgecolor','none'); axis('equal'); axis('tight');

xlabel(ax,'x [m]'); ylabel(ax,'y [m]'); title(ax,'Test pointObj, lineObj and areaObj');

HK    = gr.const(hk);
VK    = gr.const(hk);
SY    = gr.const(sy);
SS    = gr.const(ss);
STRTHD= bsxfun(@times,DEM,ones(1,1,gr.Nlay)); % starting heads

%% Boundaries of the different regions

modelArea      = area2Obj(gr,testData,'extent','NIL','name','model',DEM);
modelArea.plot(ax,'r');

outcrops       = area2Obj(gr,testData,'extent','NIL','name','outcrop');
outcrops.plot(ax,'c');

fluxBoundaries = lineObj( gr,testData,'boundaries','FLUX','type','flux');
fluxBoundaries.plot(ax,'ko');

headBoundaries = lineObj( gr,testData,'boundaries','CHD','type','head',DEM,DEM);
headBoundaries.plot(ax,'b');

genHeadBoundary =lineObj( gr,testData,'ghboundaries','GHB',DEM);
genHeadBoundary.plot(ax,'gp-.');

rivers         = lineObj( gr,testData,'rivers','RIV',DEM);
rivers.plot(ax,'b','linewidth',2);

khettaras      = lineObj( gr,testData,'khettaras','DRN','type','Qanat',DEM);
khettaras.plot(ax,'g');

wells         = pointObj( gr,testData,'wells','WEL');
wells.plot('k','markersize',8);

pumpAreas      = area2Obj(gr,testData,'areas','FLUX');
pumpAreas.plot(ax,'y');

piezom         = pointObj(gr,testData,'piezom','NIL');
piezom.fill(8);

profiles = lineObj(gr,testData,'profiles','nil');
profiles.plot(ax,'r--');


%% IBOUND array and STRTHD (starting heads)

ibound= zeros(gr.Ny,gr.Nx);           % top layer of IBOUND

ibound([modelArea.Idx]) = true;    % model area active

for ioc=1:numel(outcrops), ibound([outcrops(ioc).Idx]) = false; end

IBOUND = bsxfun(@times,ibound,ones(1,1,gr.Nlay)); % extend to all layers

%% Get recharge values from spreadsheet PER if you want to vary it in space

%% Save pass to mf_analyze
save underneath DEM testData

fprintf('***** %s -- finished. Running it took %g seconds *****\n',mfilename,toc);


