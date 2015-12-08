%% Analyzing model output
load('name'); load(basename);
load underneath  % extra parameters from mf_adapt

%% Get unformatted data
mf_checkdir;

C=readMT3D('MT3D001.UCN'); C=maskHC(C,[0,Inf]);
T=readMT3D('MT3D002.UCN'); T=maskHC(T,[0,Inf]);

%B=readBud([basename,'.bgt'],'FLOWRIGHTFACE');  % get only flow rightface
%B= mf_Psi(B);

time = [C.time]; 

%% Contour each result or make movie

crange=ContourRange(C,50);
trange=ContourRange(T,50);

figure;

% report on advection method used (and possibly also on the istherm applied)
[~,~,~,~,adv]=getExcelData(basename,'MT3D','Vertical','mixelm'); adv=ADVmethod(adv);
[~,~,~,~,iso]=getExcelData(basename,'MT3D','Vertical','isothm'); iso=isotherm(iso);

axposC = [0.1 0.10 0.8 0.35];
axposT = [0.1 0.60 0.8 0.35];

axC = axes('position',axposC,'nextplot','add','xlim',gr.xc([1 end]),'ylim',[0 C0],'xgrid','on','ygrid','on');
axT = axes('position',axposT,'nextplot','add','xlim',gr.xc([1 end]),'ylim',[0 T0],'xgrid','on','ygrid','on');

xlabel(axC,'x [m]'); ylabel(axC,'Concentration');
xlabel(axT,'x [m]'); ylabel(axC,'Temperature');

switch thisCase
    case {4,5}
    set(axC,'ylim',[Cini Cini+(C0-Cini)/20]);
    set(axT,'ylim',[Tini Tini+(T0-Tini)/20]);
end

ts1C = sprintf('MT3DMS benchmarking problem %s,Adv: %s Reaction: %s, time=%%.0f d','Conc',adv.name,iso.name);
ts1T = sprintf('MT3DMS benchmarking problem %s,Adv: %s Reaction: %s, time=%%.0f d','Temp',adv.name,iso.name);

vidObj = VideoWriter(basename);
vidObj.FrameRate = 3;
vidObj.Quality  = 80;
vidObj.open;

for it=1:length(time)

    ts2C = sprintf(ts1C,time(it));
    ts2T = sprintf(ts1T,time(it));
    
    if it==1
        htC = title(axC,ts2C);
        htT = title(axT,ts2T);
        
%        [~,hC]=contourf(axC,gr.xc,gr.zLay(:),XS(C(it).values(1,:,[1 end])),crange);
%        [~,hT]=contourf(axT,gr.xc,gr.zLay(:),XS(T(it).values(1,:,[1 end])),trange);
        hC=plot(axC,gr.xc,C(it).values(1,:,1));
        hT=plot(axT,gr.xc,T(it).values(1,:,1));
    else
        set(htC,'string',ts2C);
        set(htT,'string',ts2T);
%        set(hC,'zdata',XS(C(it).values(1,:,[1 end])));
%        set(hT,'zdata',XS(T(it).values(1,:,[1 end])));
        set(hC,'ydata',C(it).values(1,:,1));
        set(hT,'ydata',T(it).values(1,:,1));
    end
%    set(get(hC,'children'),'edgecolor','none');
%    set(get(hT,'children'),'edgecolor','none');
    
    vidObj.writeVideo(getframe(gcf));
end

vidObj.close;
%%
sigmaT = sqrt(2*DtempD/RetT*time(end));
sigmaC = sqrt(2*DmassD/RetC*time(end));

% displacement
xt = ( q - kh*dhdx )*time(end)/peff;

aT=arrow(axT, xt/RetT + [0 sigmaT], Tini+0.32*(T0-Tini)*[1 1],2,0,'color','k','linewidth',1);
aC=arrow(axC, xt/RetC + [0 sigmaC], Cini+0.32*(C0-Cini)*[1 1],2,0,'color','k','linewidth',1);

