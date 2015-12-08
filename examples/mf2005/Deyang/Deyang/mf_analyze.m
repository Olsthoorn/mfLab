%% Deyang -- artificial recharge site feb 2014
%
% TO 140226
close all

%% Load model name and the data generated and saved by mf_adapt
load('name.mat') % get basename stored in file name.mat
load(basename);  % having retrieved basename, load the data in [basename '.mat']
load underneath

fsz=10;          % set fontsize for plots

H = readDat([basename,'.hds']);
B = readBud([basename,'.bgt']);

hrange = ContourRange(H,200);
lrange = -0.5:0.05:1.5;

time = [H.time];

xlim = [441950 442400];
ylim = [3438725 3439725];

defaults = {'nextplot','add','fontsize',fsz,'xlim',xlim,'ylim',ylim};
    
%% head relative to local datum

figure('name',[sceName ' heads'],'position',screenPos(0.8));

ax=axes(defaults{:},'clim',hrange([1 end]));

xlabel(ax,'x [m]');  ylabel(ax,'y [m]');
axis(ax,'square');

vidObj = VideoWriter([sceName 'Hds']);
vidObj.FrameRate = 10;
vidObj.Quality = 80;
vidObj.open();

for it = 1:length(time)
    ttlstr = sprintf('sceName; heads, time = %.0f d',time(it));
    
    if it==1       
           ht = title(ax,ttlstr);
        [~,hh] = contourf(ax,gr.xm,gr.ym,H(it).values(:,:,end),hrange,'edgecolor','none');

        well = well.plotXY(ax,'ro');

        gePaths.plot

        hb = colorbar; set(get(hb,'title'),'string','head [m/d]');
    else
        set(ht,'string',ttlstr);
        set(hh,'zData',H(it).values(:,:,end));
        well = well.plotXY(it,'ro');
    end
    vidObj.writeVideo(getframe(gcf));
end
vidObj.close;


%% Recharge rate the leakage per m2
    
figure('name',[sceName ' infiltration'],'position',screenPos(0.8));
ax=axes(defaults{:},'clim',lrange([1 end]));

xlabel(ax,'x [m]');
ylabel(ax,'y [m]');

axis(ax,'square');

vidObj = VideoWriter([sceName 'Inf']);
vidObj.FrameRate = 10;
vidObj.Quality = 80;
vidObj.open();

for it=1:length(time)
    ttlstr = sprintf('%s, infiltration at time = %.0f d',sceName,time(it));

    GHB = B(it).term{strmatchi('HEADDEPBOUNDS',B(it).label)};
    RIV = B(it).term{strmatchi('RIVERLEAKAGE' ,B(it).label)};

    INF = sum(GHB+RIV,3)./gr.AREA;

    if it==1
           ht   = title(ax,ttlstr);
        [~,hh] = contourf(ax,gr.xm,gr.ym,INF,lrange,'edgecolor','none');
        well = well.plotXY(ax,'ro');

        gePaths.plot

        hb = colorbar; set(get(hb,'title'),'string','infiltr. [m/d]');
    else
        set(ht,'string',ttlstr);
        set(hh,'zData',INF);
        well = well.plotXY(it,'ro');
    end
    vidObj.writeVideo(getframe(gcf));
end
vidObj.close();

%% Conc

crange = 0:0.05:1;

C=[];
for iComp = numel(species):-1:1
    C{iComp} = readMT3D(sprintf('MT3D00%d.UCN',iComp));
end
well = well.setCout(C{1},C{2},C{3},B);

for iComp=1:numel(species)
    
    figure('name',[sceName ' Conc ' species{iComp} ],'position',screenPos(0.8));
    ax=axes(defaults{:},'clim',crange([1 end]));

    xlabel(ax,'x [m]');
    ylabel(ax,'y [m]');

    axis(ax,'square');

    vidObj = VideoWriter([sceName species{iComp}]);
    vidObj.FrameRate = 10;
    vidObj.Quality = 80;
    vidObj.open();

    for it=1:length(time)
        ttlstr = sprintf('%s, conc %s at time = %.0f d',sceName,species{iComp},time(it));

        c = mean(C{iComp}(it).values,3);

        if it==1
               ht   = title(ax,ttlstr);
            [~,hh] = contourf(ax,gr.xm,gr.ym,c,crange,'edgecolor','none');
            well = well.plotXY(ax,'ro');

            gePaths.plot

            hb = colorbar; set(get(hb,'title'),'string','conc [-]');
        else
            set(ht,'string',ttlstr);
            set(hh,'zData',c);
            well = well.plotXY(it,'ro');
        end
        vidObj.writeVideo(getframe(gcf));
    end
    vidObj.close();

end
%%

figure('name','wellOutput conc','pos',screenPos(0.8));
for iComp=1:numel(species)
    ax = subplot(3,1,iComp,'nextplot','add','fontSize',fsz);
    ylabel('Relative conc');
    well.plotCout(ax,species{iComp},iComp);
end
xlabel(ax,'time [d]'); 

%% Show observation wells

obsWels = pointObj(gr,basename,'Deyang3','NIL','name',{'OW','Peil'});
obsWels.plotHead(H);

%% water budget
zonebudget(B); % budget of layer 1
