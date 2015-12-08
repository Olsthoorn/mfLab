%% Analyzing model output
load('name'); load(basename); 

%% Analyzing model output

mf_checkdir();

mf_checkdir('.HDS','.UCN');

H=readDat([basename,'.hds']);
C=readMT3D('MT3D001.UCN');

hrange=   ContourRange(H,50);
crange=[0 ContourRange(C,50)];

% report on advection method used (and possibly also on the istherm applied)
[~,~,~,~,adv]=getExcelData(basename,'MT3D','Vertical','mixelm'); adv=ADVmethod(adv);
[~,~,~,~,iso]=getExcelData(basename,'MT3D','Vertical','isothm'); iso=isotherm(iso);

ts1=sprintf(['Const conc source in uniform flow field, Zheng (1999), p139, t=%%4.0f d, ',...
        'Solution method = %s'],adv.name);


%% Animation
figure; axes('nextplot','add','xgrid','on','ygrid','on','clim',crange([1 end]));
xlabel('xGr [m]');
ylabel('yGr [m]');
    
vidObj=VideoWriter(basename);
vidObj.FrameRate=3;
vidObj.open;

for it=1:length(C)
    
    ts2= sprintf(ts1,C(it).time);

    if it==1
        ht=title(ts2);
        [~,hc]=contourf(gr.xc,gr.yc,C(it).values,crange);
        colorbar;
        hdl=mf_logo; set(hdl,'color','w');
    else
        set(ht,'string',ts2);
        set(hc,'zdata',C(it).values);
    end
    set(get(hc,'children'),'edgecolor','none');
    
    vidObj.writeVideo(getframe(gcf));
end

vidObj.close;

%% Observation points

Points=[200 200; 50 50; 25 50;];         % show these three poins given their x,y
cells=([80,20,1; 60 40 1; 50 50 1]);     % show these three cells given their cell indices

% Observe what happens in above-specified points and model cells
D1=mf_observe(gr,'Points',Points,'C1, ',C,'C2, ',C,'Heads',H); % observation points given in model coordinates
D2=mf_observe(gr,'Cells' ,cells ,'C1, ',C,'Heads, ',H); % observation points given using cell indices
%D3=mf_observe(gr,'Cells' ,Cells ,{'QW','WELLS','CONSTANT HEAD'},B,'Conc Species 1',C);

figure; axes('nextplot','add','xGrid','on','yGrid','on');
xlabel('time [d]');
ylabel('Conc [mg/L]');

%% Plot the concentration development in the points
for ip=1:length(D1), plot(D1(ip).time,D1(ip).values);  end;
legend({D1.label});
title([D1.label])

%% Plot the concentration development in the cells
for ip=1:length(D2),     plot(D2(ip).time,D2(ip).values); end;
legend({'C1','C2','Head'});
ht=title([D2.label]);
set(ht,'string',['Concentration and heads on the observation points: ',get(ht,'string')]);
