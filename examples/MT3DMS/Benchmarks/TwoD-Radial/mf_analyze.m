%% Analyzing model output
load('name'); load(basename); load('underneath');  % stores the injection concetration from mf_adapt

problem='2D-radial transport';%% B=readBud([basename,'.bgt'],'FLOWRIGHTFACE');

%% Visualization
ADVmethod={'ULTIMATE (TVD)', 'FDM standard', 'MOC', 'MMOC', 'HMOC'};
isothm={'FDM standard', 'Linear', 'Freundlich', 'Langmuir','First order kinetic', 'Dual domain mass transfer', 'Dual dom. mass transf. +sorption'};

%% Analyzing model output
load('name'); load(basename);

if ~exist('H','var');  H=readDat([basename,'.hds']);   end
if ~exist('C','var');  C=readMT3D('MT3D001.UCN');      end

crange=ContourRange(C,50);

%report('mflow',basename);

[~,~,~,~,MIXELM]=getExcelData(basename,'MT3D','V','MIXELM');

%% Animation
figure; axes('nextplot','add','xgrid','on','ygrid','on','clim',crange([1 end]));
xlabel('xGr [m]');
ylabel('yGr [m]');
    
vidObj=VideoWriter(basename);
vidObj.FrameRate=3;
vidObj.open;

for it=1:length(C)
    st=sprintf(['Well in uniform flow field, Zheng (1999), p139, t=%4.0f d, ',...
        'Solution method = %s'],C(it).time,MIXELM);

    if it==1
        ht=title(st);
        [~,hc]=contourf(gr.xc,gr.yc,C(it).values,crange);
        colorbar;
        hdl=mf_logo;
    else
        set(ht,'string',st);
        set(hc,'zdata',C(it).values);
    end
    set(get(hc,'children'),'edgecolor','none');
    drawnow;
    
    vidObj.writeVideo(getframe(gcf));
end

vidObj.close;

%% Observation points

cells=[gr.Ny 1; gr.Ny 5; gr.Ny 10; gr.Ny 20];         % show these three poins given their x,y

% Observe what happens in above-specified points and model cells
D=mf_observe(gr,'Cells' ,cells ,'C, ',C); % observation points given using cell indices

figure; axes('nextplot','add','xGrid','on','yGrid','on');
xlabel('time [d]');
ylabel('Conc [mg/L]');

for ip=1:length(D), plot(D(ip).time,D(ip).values);  end;
legend({D.label});

%% Graph along an axis
figure, 
axes('nextplot','add','xgrid','on','ygrid','on');
title([problem,sprintf(' t=%g d, TVD',C(end).time)]);
xlabel('radial distance [m]'); ylabel('C/Co');

plot(gr.xc,C(end).values);
