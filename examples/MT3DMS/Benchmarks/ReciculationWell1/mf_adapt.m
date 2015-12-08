%% Simulating the effect of a recirculation wells (Zheng, 2005, MT3DMS v5.3, p33..)
%
% The recirculation well works, and also the multi-node well MNW1 works (by
% switching off WEL and switchting on MNW1 in the NAM worksheet). But the
% recirculation well does not work with the multi-node well, as becomes
% evident from the ouput. It does not matter whether the multi-node well
% has more than one cell per screen, as is the case in the example
% RecirculationWell2, in this exampleRecirculationWell1, each well has only
% one cell, and it still does not work. This is an issue in the Fortran
% code of the packages and its support by Mt3DMS.
%
% TO 121126 121128

clear variables; close all;

BACKGROUND=true;

basename='RecirculationWell1';

[hdr,data]=getExcelData(basename,'Data','vert');

iswell = logical(data(strmatchi('iswell',hdr)));

%% Grid
xGr = 100*(0:1:46);
yGr = 100*(0:1:31);
zGr = [0 -50];

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
if iswell
    well=wellObj(basename,'wells',gr,HK,{'PER' 'Q'},{'PER','C'});
else
    well=MNW1Obj(basename,'wells',gr,HK,{'PER' 'Q'},{'PER','C'});
end

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

save underneath dHdx dHdy iswell
