% MF_ADAPT: Example a rainwater lens in a shallow Holocene cover layer with seepage from a regional aquifer
%
% USAGE:
%   mf_adapt  % to juse see of it works
%   mf_setup  % in invoke it and run the models
%             % when model has finihsed run
%   mf_analyze
%
%% Build-up of a rainwater lense caused by precipitaton surplus, on average.
% The cross section is flanked by two ditches with fixed heads. At the
% bottom of the cross section the head is maintained at the level of the
% two ditches using CHD. Rain has concentration 0, the lower boundary has
% concentration cGr.
%
% TO 110511 120514

clear variables; close all;

basename='RainLense';

GREP='STRESS';

%% The Model

Evapotranpiration_factor = 1;    % to allow playing arouond with relative
                                 % evaportranspiataion as opposed to precipitation.
                                 % Set 1 for default

cGrw =1;
cfresh=0;
peff=0.35;

iDitch = 2; % Ditch as zone 
iCHD   = 3; % CHD boundary as zone

%% Grid
xGr=-50:1:50; 
yGr=[0 1];
zGr=[0 -5:-5:-50];  % Plane elevation vector

gr=gridObj(xGr,yGr,zGr);

%% Basic arrays
ditchbot=-1;

[HK              ]=mf_zone(basename,xGr,yGr,zGr,'Config','Material','kh');
[VK              ]=mf_zone(basename,xGr,yGr,zGr,'Config','Material','kv');
[PEFF            ]=mf_zone(basename,xGr,yGr,zGr,'Config','Material','peff');
[STCONC          ]=mf_zone(basename,xGr,yGr,zGr,'Config','Material','stconc');
[STRTHD          ]=mf_zone(basename,xGr,yGr,zGr,'Config','Material','strthd');
[SS              ]=mf_zone(basename,xGr,yGr,zGr,'Config','Material','ss');
[SY,Conf,material]=mf_zone(basename,xGr,yGr,zGr,'Config','Material','sy');

STRTHD(:,:,:)=0;

IBOUND=gr.const(1); IBOUND(:,[1 end],1)=-iDitch;
IBOUND(:,:,end) = -iCHD;

ICBUND=IBOUND;

%% Speciy infiltration n basis as recharge within barrier

[~,~,NPER]=getPeriods(basename);

%% Read precipiptation - evapotranspiration file

fid=fopen('PE-92-00.txt','r');

pne=fscanf(fid,'%d-%d-%d %f %f',[5,Inf])';
pne=[datenum(pne(:,3),pne(:,2),pne(:,1)),pne(:,end-1:end)/1000];

fclose(fid);

%% Change to weekly data to redue amount of data

RECH=pne(:,2)-Evapotranpiration_factor*pne(:,3); % <--- recharge

nweek=fix(length(RECH(:))/7);

RECH=mean(reshape(RECH(1:nweek*7),[nweek,7]),2);

[CHD,PNTSRC] = bcnZone(basename,'CHD',abs(IBOUND),{iCHD STRTHD(IBOUND==-iCHD) STRTHD(IBOUND==-iCHD)},cGrw);
    
%% set the wells

%[wel,WEL,PNTSRC]=gr.well(basename,HK,'wells');

save underneath Conf
