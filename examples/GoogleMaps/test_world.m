%% World figure at zoom level 0, with locations plotted

locations={...
[ 25.253159, 55.178000], 'Dubai-world island';
[ 52.372821,  4.893667], 'Dam Monument Amsterdam';
[ 68.962970, 33.089563], 'Murmansk, Russia';
[-33.856857,151.215192], 'Sydney';
[-54.795444,-68.232218], 'Ushuaia, Argentina';
[ 37.808810,-122.409803],'San Francisco, Fisherman''s Wharf';
[ 40.748524,-73.985676], 'New York Empire State building'
};

%% GM object
center = [0 0];
zoom   = 0;
pixels = [256;256];
maptype= 'satellite';

GM = googleMapObj(center,zoom,pixels,maptype,'png');
GM.image([0 1],[1 0]);
set(gca,'ydir','reverse');

%%
% reposition and resize figure
scr=get(0  ,'screensize'); sw= scr(3); sh=scr(4);
set(gcf,'position',[scr(3)/6,scr(4)/6,2/3*scr(3),2/3*scr(4)]);

%% Locations
locs = cat(1,locations{:,1});
gp = googlePointObj(locs(:,1),locs(:,2),zoom);
plot([gp.x],[gp.y],'o','markerfacecolor','y');

offset=[0.01,-0.02];
for iLoc=1:size(locations,1)
    text(gp(iLoc).x+offset(1),gp(iLoc).y+offset(2),locations{iLoc,end},'color','w');
end

%% Lets splot the merdians and latitudes

%%
% Meridians and latitudes
N= [85 80:-10:-80 -85];
E=-180: 20:180;

%%
% in mercator x and y
yN = mf_GMlat2y(N);
xE = mf_GMlon2x(E);

%%
% set tick of axes equal to meridians and latitudes
set(gca,'yTick',yN,'yTickLabel',N,'xgrid','on','ycolor',[ 1 1 1]);
set(gca,'xTick',xE,'xTickLabel',E,'ygrid','on','xcolor',[ 1 1 1]);

%%
% background color of figure
set(gcf','color',[0.3 0.3 0.3]);

title('Locations in the previous figure','color','y');
