Seepage = 0.0005;

N = TPE(:,2)-TPE(:,3);
t = TPE(:,1); dv = datevec(t); winter=dv(:,2)>=10 & dv(:,2)<=3;
dt= diff(t);
dt = [dt(1); dt];

sigmaMu = 0.40;
muMax   = 0.10;
h50     = hmv-0.35; 
mu = @(h) 0.5 * muMax  * erfc((h - h50)/sigmaMu);
hGreppels = -4.45; % Afvoer via het taludvlak
hBoezem   = -4.80;
Twin=  7.0;
Tsun=  100;
hmv = -3.95;
h = zeros(size(t));
for it=1:numel(t)
    if it==1, hPrev = hsl; else hPrev = h(it-1); end
    if hPrev>hGreppels
        hSl = hGreppels;
        T=Twin;
        gamma = T/mu(hPrev);
    else
        hSl = hBoezem;
        T=Tsun;
        gamma=T/mu(hPrev);
    end
    h(it) = hSl + (hPrev-hSl)*exp(-dt(it)/T) + (N(it)+Seepage)*gamma .* (1-exp(-dt(it)/T));
    fprintf('%3d %s %12.3f %12.4f\n',it,datestr(t(it)),h(it),N(it));
    h(it) = min(hmv,h(it));
end

if true
plot(t(I),h(I),'m');
else
figure; 
tStart = datenum(2010,1,1); tEnd=datenum(2010,12,31);
axes('nextPlot','add','xGrid','on','yGrid','on','color','none','fontsize',12,...
    'xlim',[tStart tEnd],'ylim',[-6 -3.5]);
xlabel('tijd');
ylabel('hNAP');
I=t>=datenum(2010,1,1);
plot([tStart tEnd],[hmv hmv],'g');
plot(t(I),h(I),'k','linewidth',2);
datetick;
end