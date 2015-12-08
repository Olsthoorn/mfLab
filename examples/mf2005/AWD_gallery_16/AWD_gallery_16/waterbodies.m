% Prepare the waterbodies for this model
% done in separate m-script to reduce complexity of body of model
% just start with load Waterbody after having successfully run this file
%
%
% TO 100922

load('model_Muriel/DR16.mat');  % loading rch_pond drain canal VBAK

Waterbody=[];
k=0;
for i=1:length(KA)
    k=k+1;
    Waterbody(k).name=KA(i).NAAM;
    Waterbody(k).shortnm=KA(i).NAAM;
    Waterbody(k).type='Canal';
    Waterbody(k).code=KA(i).CODE;
    Waterbody(k).nr=KA(i).KANAAL_AWD;
    Waterbody(k).gauge=KA(i).PEILSCHAAL;
    Waterbody(k).x=KA(i).x;
    Waterbody(k).y=KA(i).y;
    Waterbody(k).area=KA(i).AREA;
    Waterbody(k).perimeter=KA(i).PERIMETER;
    Waterbody(k).start=datenum(KA(i).START  ,1,1);
    Waterbody(k).end  =datenum(KA(i).EINDE+1,1,1);
    Waterbody(k).h_avg=NaN;
    Waterbody(k).cbot=NaN;
end

for i=1:length(GE)
    k=k+1;
    Waterbody(k).name=GE(i).NAAM;
    Waterbody(k).shortnm=Waterbody(k).name;
    Waterbody(k).type='Pond';
    Waterbody(k).code=GE(i).CODE;
    Waterbody(k).nr  =GE(i).GEULEN_AWD;
    Waterbody(k).gauge=NaN;
    Waterbody(k).x=GE(i).x;
    Waterbody(k).y=GE(i).y;
    Waterbody(k).area=GE(i).AREA;
    Waterbody(k).perimeter=GE(i).PERIMETER;
    Waterbody(k).start=datenum(KA(i).START  ,1,1);
    Waterbody(k).end  =datenum(KA(i).EINDE+1,1,1);
    Waterbody(k).h_avg=GE(i).GEM_PEIL_7;
    Waterbody(k).cbot=GE(i).C_BODEM;
end

spg_pond.name='spg_pond';
[spg_pond.x,spg_pond.y]=kmlpath('DR16SpgPond.kml'); % Hand digitized in Google Earth

for i=1:length(spg_pond)
    k=k+1;
    Waterbody(k).name=spg_pond.name;
    Waterbody(k).shortnm=Waterbody(k).name;
    Waterbody(k).type='spg_pond';
    Waterbody(k).code=NaN;
    Waterbody(k).nr  =NaN;
    Waterbody(k).gauge=NaN;
    Waterbody(k).x=spg_pond.x;
    Waterbody(k).y=spg_pond.y;
    Waterbody(k).area=NaN;
    Waterbody(k).perimeter=NaN;
    Waterbody(k).start=datenum(1965,1,1);
    Waterbody(k).end  =datenum(9999,1,1);
    Waterbody(k).h_avg=NaN;
    Waterbody(k).cbot=NaN;
end

for i=1:length(Waterbody)
    Waterbody(i).xCtr=mean(Waterbody(i).x);
    Waterbody(i).yCtr=mean(Waterbody(i).y);
end

shortnames={'Hup_1','Hup_2','Hup_3','Hup_4','PK15','G27',...
'G28N','G28S','Toev_1','Toev_2','spg_pond'};  % this is hard-wiring and therefore order dependent!!

for i=1:length(Waterbody)
    Waterbody(i).shortnm=shortnames{i};
end


%% Cleanup the drain

% only need center of VBAK (central shaft of drain)
shaft.x=mean(VBAK(:,1));
shaft.y=mean(VBAK(:,2));
shaft.z=mean(VBAK(:,3));

% remove all intermediate points of drain in shape, only keep its ends
% just merging and retaining min and max is ok because drain is orientetad SW NE
xDr=[min([drain.x]) max([drain.x])];
yDr=[min([drain.y]) max([drain.y])];

% first shift the connection shaft to drain
[x_shaft,y_shaft]=point2line(xDr,yDr,shaft.x,shaft.y);  % dr_shaft comes from load DR16_input.mat
 z_shaft=shaft.z;

 dr(1).name='DR16south';
 dr(2).name='DR16north';

dr(1).x=[xDr(1) x_shaft];
dr(1).y=[yDr(1) y_shaft];
dr(1).z=zeros(size(dr(1).x));
dr(1).xshaft=x_shaft;
dr(1).yshaft=y_shaft;
dr(1).zshaft=z_shaft;
dr(1).rw=0.15;       % drain radius
dr(1).Rw=0.45;       % gravel pack radius


dr(2).x=[x_shaft xDr(end)];
dr(2).y=[y_shaft yDr(end)];
dr(2).z=zeros(size(dr(2).x));
dr(2).xshaft=x_shaft;
dr(2).yshaft=y_shaft;
dr(2).zshaft=z_shaft;
dr(2).rw=0.15;       % drain radius
dr(2).Rw=0.45;       % gravel pack radius

drain=dr;

save Waterbodies Waterbody drain


