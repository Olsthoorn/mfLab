%% Analyzing model output
load('name'); load(basename); load('underneath');

C{1}= readMT3D('MT3D001.UCN');
C{2}=readMT3D('MT3D001S.UCN');

%[~,~,~,~,NCOMP]=getExcelData(basename,'MT3D','V',NCOMP);

%% Visualization

Points=[8, 0.5];

% report on advection method used (and possibly also on the istherm applied)
[~,~,~,~,adv]=getExcelData(basename,'MT3D','Vertical','mixelm'); adv=ADVmethod(adv);
[~,~,~,~,iso]=getExcelData(basename,'MT3D','Vertical','isothm'); iso=isotherm(iso);

ts1=sprintf(['Const conc source in uniform flow field, Zheng (1999), p139, t=%%4.0f d, ',...
        'Solution method = %s'],adv.name);

Obs=mf_observe(gr,'points',Points,'Dissolved, ',C{1},'Sorbed, ',C{2});

st={'Concentration dissolved','Concentration sorbed'};
sl={'Dissolved conc [mg/L', 'Sorbed cocc [mg/g]'};

figure;

for iObs=1:length(Obs)
    subplot(length(Obs),1,iObs,'nextplot','add','xgrid','on','ygrid','on');
    xlabel('time [min]');
    ylabel(sl{iObs});
            
    s=sprintf('Mt3D benchmark: %s, %s, x=%d cm, advmethod=%s, isotherm=%s',...
        basename,sl{iObs},Points(1),adv.name,iso.name);
        
    title(s);

    plot(Obs(iObs).time,Obs(iObs).values);
end

