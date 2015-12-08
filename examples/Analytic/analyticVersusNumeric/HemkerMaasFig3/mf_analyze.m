%% Analyzing output of the model
% TO 091011 091129 120413
 
load('name.mat') % get basename stored in file name.mat
load(basename);  % having retrieved baasename load the data in basename.mat
load underneath  % to get gr object

defaults =  {'nextplot','add','ydir','reverse','xgrid','on','ygrid','on'};

IL = find(HK(1,1,:)>1);

drawdown=hantushn(Q,gr.xm,t,St,Sf(2:end-1),c,T(2:end-1));

H=readDat([basename,'','.hds']); % read the unformatted head file

figure('name','Hemker Maas (1987) fig3a,b','pos',screenPos(0.75));
ax(1) = subplot(1,2,1,defaults{:},'xscale','lin','yscale','lin','xlim',[0 6000],'ylim',[0 1]);
ax(2) = subplot(1,2,2,defaults{:},'xscale','log','yscale','log','xlim',[10 1e4],'ylim',[1e-3 10]);
title(ax(1),'fig HM 3a');
title(ax(2),'fig HM 3b');
xlabel(ax(1),' r --> m');
xlabel(ax(2),' r --> m');
ylabel(ax(1),' s [m]');
ylabel(ax(2),' s [m]');

clrs = 'brgmk';
leg=[]; % Nlay = size(drawdown,1);
Nlay=2;
for it=1:8:numel(t)
    for iL=Nlay:-1:1
        plot(ax(1),gr.xm,XS(H(it).values(:,:,IL(iL))),['o' mf_color(it,clrs) '--']);
        plot(ax(2),gr.xm,XS(H(it).values(:,:,IL(iL))),['o' mf_color(it,clrs) '--']);
        plot(ax(1),gr.xm,drawdown(iL,:,it),[mf_color(it,clrs) '-'],'lineWidth',2);
        plot(ax(2),gr.xm,drawdown(iL,:,it),[mf_color(it,clrs) '-'],'lineWidth',2);
        leg{(it-1)*Nlay+iL} = sprintf('%e/%d',t(it),iL);
    end
end
%legend(ax(1),leg);
%legend(ax(2),leg);

%% Plot the dradown versus time for the piezometer at r

r = 10; ir=hit(gr.xGr,r);

drawdown=hantushn(-Q,r,t,St,Sf(2:end-1),c,T(2:end-1));
dd   = abs(squeeze(drawdown));

ddMF = NaN(size(dd));
for it=1:numel(H), ddMF(:,it) = XS(H(it).values(1,ir,IL)); end

ttl = sprintf('HM (1987) problem fig 3, %d sublayes per aquitard',layersPerAquitard);
ddMF = abs(ddMF);

figure('name','HM fig3, vs time');

ax(3) = subplot (2,1,1,defaults{:},'xscale','log','yscale','lin','ylim',[1e-4 10]);
ax(4) = subplot (2,1,2,defaults{:},'xscale','log','yscale','log','ylim',[1e-4 10]);

xlabel(ax(3),'time [d]'); ylabel(ax(3),'dd[m]');
xlabel(ax(4),'time [d]'); ylabel(ax(4),'dd[m]');
title(ax(3),ttl);         title(ax(4),ttl);
plot(ax(3),t,dd,'r','linewidth',2);
plot(ax(3),t,ddMF,'bo-');
plot(ax(4),t,dd,'r','linewidth',2);
plot(ax(4),t,ddMF,'bo-');

plot(t,ddMF,'o-b');
plot(ax(4),t,ddMF,'o-b');
