%% Simulating the effect of a recirculation wells (Zheng, 2005, MT3DMS v5.3, p33..)
%
% Same example as RecirculationWell1, but more cells (layers) per well,
% real 3D.
%
% TO 121126

clear variables; close all;

BACKGROUND=1;

basename='RecirculationWell2';

%% Grid
xGr = 100*(0:1:46);
yGr = 100*(0:1:31);
zGr = 0:-5:-150;

gr=gridObj(xGr,yGr,zGr);

%% 3D arrays
HK    = gr.const(10);
VK    = gr.const(10);
PEFF  = gr.const(0.35);
IBOUND= gr.const(1); IBOUND([1 end],:,:)=-1; IBOUND(:,[1 end],:)=-1;
STCONC= gr.const(0);

% Starting heads with gradient
dHdx = -0.003; dHdy = 0; strthd =  0;
STRTHD = strthd+(gr.XM-mean(gr.xGr))*dHdx+(gr.YM-mean(gr.yGr))*dHdy;

% WELLS
well=MNW1Obj(basename,'wells',gr,HK,{'PER' 'Q'},{'PER','C'});
%well=wellObj(basename,'wells',gr,HK,{'PER' 'Q'},{'PER','C'});

well=well.screenPoint(gr,0.5);

NPER=numel(well(1).Q);

% put recirculation wells
for iw=numel(well):-1:1
    iPeer = well(iw).UserData.recircWellNr;
    peer = well([well.nr]==iPeer);
    if ~isempty(peer)
        well(iw).C = -ones(1,NPER)*peer.UserData.screenPoint.idxMT3D;
    end
end

ICBUND = IBOUND;

save underneath dHdx dHdy
