%% Analyzing model output
load('name'); load(basename);

% pack C into cell to deal with more than one component or both dissolved
% and sorbed concentrations
C{1}= readMT3D('MT3D001.UCN');
C{2}=readMT3D('MT3D001S.UCN');

%[~,~,~,~,NCOMP]=getExcelData(basename,'MT3D','V',NCOMP);

%% Visualization

Points=[8, 0.5];

% report on advection method used (and possibly also on the istherm applied)
[MT3nams,MT3vals,MT3thdr,MT3txt,adv]=getExcelData(basename,'MT3D','Vertical','mixelm');
[~,~,~,~,iso]                       =getExcelData(basename,'MT3D','Vertical','isothm');

adv=ADVmethod(adv);
iso=isotherm(iso);

Obs=mf_observe(gr,'points',Points,'Dissolved, ',C{1},'Sorbed, ',C{2});

st={'Concentration dissolved','Concentration sorbed'};
sl={'Dissolved conc [mg/L', 'Sorbed cocc [mg/g]'};

figure;

for iObs=1:length(Obs)
    subplot(length(Obs),1,iObs,'nextplot','add','xgrid','on','ygrid','on');
    xlabel('time [min]');
    ylabel(sl{iObs});
            
    s=sprintf('Mt3D benchmark: %s, %s, x=%d cm, advmethod=%s, isotherm=%s,',...
        basename,sl{iObs},Points(1),adv.name,iso.name);
        
    title(s);

    plot(Obs(iObs).time,Obs(iObs).values);
end

