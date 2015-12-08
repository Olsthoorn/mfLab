%% Select PTE and heads

%readWarecoPB --> pbuis(1).TH, pbuis(2).TH, pbuis(3).TH
% PB10 = pbuis(1).TH;
% PB11 = pbuis(2).TH;
% PB12 = pbuis(3).TH;

tStart= datenum(2010,4,10);
tEnd  = datenum(2010,4,20);
tStart= datenum(2010,8,1);
tEnd  = datenum(2010,8,31);

I = TPE(:,1)>=tStart & TPE(:,1)<=tEnd;

J = PB10(:,1)>=tStart & PB10(:,1)<=tEnd;
K = PB11(:,1)>=tStart & PB11(:,1)<=tEnd;
L = PB12(:,1)>=tStart & PB12(:,1)<=tEnd;

defaults = {'nextPlot','add','xGrid','on','yGrid','on','xLim',[tStart tEnd]};

figure;
ax1 = axes(defaults{:}); datetick;
xlabel(ax1,'2010');
ylabel(ax1,'P & E m/d');
ax2 = axes('position',get(ax1,'position'),defaults{:},'color','none','yAxisLocation','right');
set(ax1,'xlim',[tStart tEnd]); datetick(ax1);
set(ax2,'xlim',[tStart tEnd]); datetick(ax2);
datetick(ax2)
ylabel(ax2,'m NAP');

bar(ax1,TPE(I,1),TPE(I,2),'b','lineWidth',2);
bar(ax1,TPE(I,1),TPE(I,3),'r','lineWidth',2);
plot(ax2,PB10(J,1),PB10(J,2),'g');
plot(ax2,PB11(J,1),PB11(J,2),'b');
plot(ax2,PB12(J,1),PB12(J,2),'r');