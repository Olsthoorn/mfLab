%% Analyzing output of xsecton model
% TO 110423
 
clear variables
close all

%% load model name and basename contained in name.mat
load name
load(basename);  % this yields the stored matrices that make up the grid and the model
load underneath; % additional values from mf_adapt

startyear=0; % start year

aday=24; % start year and days per year

ttl='cross section ';

%% load the unformatted files with the heads, the concentration and the% budget terms

H=readDat([basename '.hds']); H=maskHC(H,IBOUND);   % read heads and mask inactive cells with NaN
C=readMT3D('MT3D001.UCN');    C=maskHC(C,ICBUND);   % read concentrations and mask inactive with NaN

B=readBud([basename '.bgt']);  B=maskHC(B,IBOUND); % read budgetfile and mask
B=mf_Psi(B); % fpritnf('Stream function added to budget struct

%%

extr=zeros(size(B));
for i=1:length(B)
    extr(i)=sum(B(i).term{strmatchi('WELL',B(i).label)}(:));
end

%% make movie

dc  =0.01; crange  =contourrange(C,dc);
dpsi=0.5; psirange=contourrange(B,dpsi,'','Psi');
dphi=0.5; phirange=contourrange(H,dphi);
dvk= 0.1; vkrange =contourrange(log(VK),dvk);

%% Heade figure
figure('position',screenPos(0.75));

bgr=get(gcf,'color');
headfig.ax(1)=axes('position',[0.1 0.81 0.8 0.1],'color','none');
headfig.ax(2)=axes('position',[0.1 0.81 0.8 0.1],'color','none');

headfig.A{1}=imread('Ridderkerk1.png'); % replace by better picture showing river and surrounding polder
headfig.A{2}=imread('Ridderkerk1.png'); % replace by better picture

axes(headfig.ax(1)); headfig.hdl(1)=image(xGr,[6 -6],headfig.A{1});
axes(headfig.ax(2)); headfig.hdl(2)=image(xGr,[6 -6],headfig.A{2});

title('Schematic of subsurface storage of pumpingstation Reijerwaard, Ridderkerk, NL');
set(headfig.ax(1),'ydir','normal','xcolor',bgr,'ycolor',bgr,'color','none');
set(headfig.ax(2),'ydir','normal','xcolor',bgr,'ycolor',bgr,'color','none');

%% Other axies

pos= [0.1 0.1 0.8 0.65];

ax=[axes('position',pos);...
    axes('position',pos);...
    axes('position',pos);...
    axes('position',pos)];

cL=[ crange(  [1 end]);...
     psirange([1 end]);...
     phirange([1 end]);...
     vkrange( [1 end])
   ];

%mf_setmulticolormap(ax,cL,64);

ranges={crange,psirange,phirange}; submaplengths=64; mapnames='jet';

cmap=mf_setColormap(ax,{crange,phirange,psirange,vkrange},64,'jet');

%% Switch to new axis and print filled and unfilled contours

xlabel(ax(1),'Distance [m]'); ylabel(ax(1),'elevation [m]');
ht=title([ttl sprintf(', day %.2f  dphi=%.2f m, dPsi=%.2f m2d/d',H(i).totim/aday+startyear,dphi,dpsi)]);

[~,hc1]=contourf(ax(1), xm, squeeze(zm), XS(C(1).values) ,crange,'edgeColor','none');  % ,'linecolor','none');

% axes(ax(2));  [c,hp] =contour( xm, zm(:), XS(H(i).values), phirange,'b');
[~,hs] =contour(ax(3), xGr(2:end-1), zGr(:)', B(1).Psi, psirange,'color',[0.4 0.4 0.4]);

%% Plot the aquitards, transparently
plotConf(Conf,{'K1','K2'},[0.2 0.5 0.2],0.3);

linkaxes(ax);

set(ax(1),'xlim',[xGr(1) xGr(end)]);

%% plot wells but do not show yet

d=3;
for iw=1:length(well)
       plot(well(iw).x+[-d -d d d],[zGr(1) well(iw).z([2 2]) zGr(1)],'color',bgr');
       well(iw).hdl=fill(well(iw).x+[-d d d -d],well(iw).z([1 1 2 2]),'w');
       set(well(iw).hdl,'visible','off');
end

%% loop through the figures to make movie

aviobj=avifile(basename,'compression','none',...
    'fps',5,'colormap',colormap,...
    'videoname','Oasen');

for i=1:length(H)

   set(hc1,'zdata',XS(C(i).values)); 
   set(get(hc1,'children'),'edgecolor','none');

%  set(hp, 'zdata',XS(H(i).values));
   set(hs, 'zdata',B(i).Psi);
   set(ht, 'string',...
       [ttl sprintf(', day =%.2f, dphi=%.2f m, dPsi=%.2f m2d/d',...
       H(i).totim,dphi,dpsi)]);
      
   drawnow;
 
    if film, F=getframe(gcf); aviobj=addframe(aviobj,F); end
end
    
aviobj=close(aviobj);

if ismac, avi_compress(basename); end


%% Extracted concentration

c=NaN(length(B),length(well));
qq=NaN(size(well));
cc=NaN(size(well));

for i=1:length(B)
    for iw=1:length(well)
        RCL=cellIndices(well(iw).idx,size(IBOUND),'RCL');
        well(iw).idxL=cellIndex([RCL(:,2)-1 RCL(:,1) RCL(:,3)],size(IBOUND));
        well(iw).idxR=cellIndex([RCL(:,2)+1 RCL(:,1) RCL(:,3)],size(IBOUND));
               
        q_=B(i).term{strmatchi('WELLS',B(i).label)}(well(iw).idx)';
        c_=0.5*(C(i).values(well(iw).idxL) + C(i).values(well(iw).idxR));
        qq(iw)=sum(q_);
        c(i,iw)=sum(c_.*q_)./sum(q_);
    end
end

figure; plot([C.time]'/365.25,c); xlabel('time [yrs]'); ylabel('Cl^- concentration'); grid on
title('Concentration in well');
legend({well.name},3);
