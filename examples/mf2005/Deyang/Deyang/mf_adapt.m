%
basename = 'Deyang';

save name basename

GREP = 'STRESS PERIOD';

scenario = 3;
switch scenario
    case 1, sceName ='Deyang1';
    case 2, sceName ='Deyang2';
    case 3, sceName ='Deyang3';
    otherwise
        error('unknown case');
end

GE      = '../GEarth/';
gePaths = kmlPathsObj([GE sceName '.kml']);
gePaths(strmatchi('Wells',{gePaths.name}))=[];

SiteZero     = gePaths(strmatchi('SiteZero',{gePaths.name}));
DeyangSite   = gePaths(strmatchi('DeyangSite',{gePaths.name}));

xMid = (min(DeyangSite.X)+max(DeyangSite.X))/2;
yMid = (min(DeyangSite.Y)+max(DeyangSite.Y))/2;
xMin  = SiteZero.X-1500;
xMax  = SiteZero.X +450;
yMin  = SiteZero.Y-1500;
yMax  = SiteZero.Y+1500;
xMinSite = min(DeyangSite.X);
xMaxSite = max(DeyangSite.X);
yMinSite = min(DeyangSite.Y);
yMaxSite = max(DeyangSite.Y);

% smooth grid using normalSpace
dxRatio= 20; sigmaX = 500;
dyRatio= 20; sigmaY = 500;

xGr = [normalSpace(xMin,xMinSite,xMinSite,25,dxRatio,sigmaX) ...
       xMinSite:5:xMaxSite ...
       normalSpace(xMaxSite,xMaxSite,xMax,20,dxRatio,sigmaX)];
yGr = [normalSpace(yMin,yMinSite,yMinSite,25,dyRatio,sigmaY) ...
       yMinSite:5:yMaxSite ...
       normalSpace(yMaxSite,yMaxSite,yMax,25,dyRatio,sigmaY)];

zSite    =     + 478;
zRoad    =  4  + zSite; % Datum
zDike    =  6  + zSite;
zWater   = -2  + zSite;
zRivBot  = -6  + zSite;
zSite    =  0  + zSite;
zBotAquif= -26 + zSite;
zWetland = +0  + zSite;

peff = 0.5;
kD   = 250;
kh   = kD/(zRoad - zBotAquif);
cRiver = 100; % d
cBasin = 2.5;  % d

xLim   = [-350 1000];
yLim   = [-1000 1000];

%zGr = [zRoad zDike zWater zRivBot zSite, zBotAquif zWetland -10:2:10];
zGr = [zDike 472 466.8 zBotAquif];
gr = gridObj(xGr,yGr,zGr);

% Arrays
IBOUND = gr.const(1);
HK     = gr.const(kh);
VK     = HK;
STRTHD = gr.const(0);

species = {'RivWater' 'RchWater' 'IrrWater'};

% Line and area2Obj
MiyungA = area2Obj(gr,basename,sceName,'RIV','name','MiyungA','conc',species);
MiyungB = area2Obj(gr,basename,sceName,'RIV','name','MiyungB','conc',species);
Basin   = area2Obj(gr,basename,sceName,'GHB','name','Basin',  'conc',species);
Creek   = lineObj (gr,basename,sceName,'RIV','name','Creek',  'conc',species);

well = wellObj(basename,'wells',gr,HK,'PER');
well = well(strcmpi(sceName,{well.name}));
% replace scenario in well.name to 'well+nr'
for iw=1:numel(well)
    well(iw).name= sprintf('well%d',iw);
end

ICBUND = IBOUND;
PEFF   = gr.const(peff);
STCONC{1} = gr.const(0);
STCONC{2} = gr.const(0);
STCONC{3} = gr.const(0);

%% Modpath

save underneath gePaths sceName zSite zDike zBotAquif species