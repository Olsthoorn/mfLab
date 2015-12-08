% drain = lineObj(basename,'khettaras',gr,'DRN',DEM);
% riv   = lineObj(basename,'river'    ,gr,'RIV',DEM,[],DEM);
% chd   = lineObj(basename,'headBoundary',gr,'CHD',DEM,DEM);
% flux  = lineObj(basename,'fluxBoundary',gt,'FLUX');
% ghb   = lineObj(basename,'ghbBounday', gr,'GHB',DEM);

minDZ=0.1;
xGr = 0:1000;
yGr = 1000:-100:0;
zGr = [0 -25 -50 -75];

gr = gridObj(xGr,yGr,zGr,'minDZ',minDZ);

%figure('name','testLineObj','pos',screenPos(0.75));
%ax = axes('nextplot','add','xlim',gr.xGr([1 end]),'ylim',gr.yGr([end 1]));
[xx yy] = ginput();

fprintf('%g\t%g\n',[xx yy]');

