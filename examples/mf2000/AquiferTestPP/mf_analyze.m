%% Analyzing model output
load('name');
load(basename);
load('underneath');

%% Known transmissivity
kDknown = sum((IBOUND(:,end,:)~=0).*HK(:,end,:).*gr.DZ(:,end,:),3);
Ssknown = sum((IBOUND(:,end,:)~=0).*SS(:,end,:).*gr.DZ(:,end,:),3);

% characteristic constant for elastic storage and specific yield respectively
TE = gr.xGr(2).^2*sum((IBOUND(:,end,2:end)~=0).*SS(:,end,2:end).*gr.DZ(:,end,2:end),3)./(4*kDknown);
TS = gr.xGr(2).^2*sum((IBOUND(:,end,1    )~=0).*SS(:,end,1    ).*gr.DZ(:,end,1    ),3)./(4*kDknown);

%% get drawdowns
H = maskHC(readDat([basename '.DDN']),[-Inf 1e4],[NaN NaN]);

% Compute drawdowns as average over the screen length in the well bore
DDN = NaN(size(kDknown,1),length(H));
for it=1:length(H)
    DDN(:,it) = mean(H(it).values(:,1,well(1).iLay),3);
end

%% Plot the drawdown
figpos = get(0,'screensize').*[1 1 0.75 0.75];
figure('position',figpos); hold on; xlabel('time [d]','fontsize',12); ylabel('drawdown [m]','fontsize',12);
title(sprintf('Computed drawdown all scenario cases, Q=%.0f m3/d',well(1).Q),'fontsize',12);
plot([H.totim]',DDN'); set(gca,'xscale','log','xgrid','on','ygrid','on','xMinorGrid','on','fontsize',12);

%% Compute drawdown increase per log cycle
time = [H.totim];

DE = NaN(size(TE),2);  % elastic  delta DDN (early time drawdown curve)
DS = NaN(size(TS),2);  % phreatic delta DDN (late time drawdown curve) 
% for i=1:length(kDknown)
%     DE(i,:) = interp1(time,DDN(i,:),[max(time(  1)   ,TE(i)*100),max(time(  1),TE(i)*1000)]);
%     DS(i,:) = interp1(time,DDN(i,:),[min(time(end)/10,TS(i)*100),min(time(end),TS(i)*1000)]);
% end

for i=1:length(kDknown)
    DE(i,:) = interp1(time,DDN(i,:),[1e-6 1e-5]);
    DS(i,:) = interp1(time,DDN(i,:),[1e+2 1e+3]);
end

%%
kDE = -log(10) * well(1).Q./(4*pi*diff(DE,1,2));
kDS = -log(10) * well(1).Q./(4*pi*diff(DS,1,2));


%% Drawdown with distance for the last time

figure; hold on; xlabel('r [m]'); ylabel('drawdown [m]'); title('Drawdown versus distance');
for iy=1:gr.Ny
    plot(gr.xm,XS(H(end).values(iy,:,gr.zm<-49 & gr.zm>-51)),mf_color(iy,'rbgkcm'));
end
grid on;
legend('1','3','10','30','100',2);
set(gca,'xscale','log','ydir','reverse');

%% Point extraction
% r= abs(diff(well(1).z)/2)/3.5;
% beta=sqrt(ss(1)/k);
% k = HK(1);
% 
% s = -well(1).Q/(2*pi*k)*erfc(beta*r/2./sqrt(time))/r;
% 
% plot(time,s,'k','linestyle','-','linewidth',2);


%% Theis without partial penetration
rw = gr.xGr(2);

s = -well(1).Q./((4*pi*kDknown)*ones(size(time))).*expint(rw^2*(Ssknown./kDknown)*(1./time)/4);

%figure;
plot(time,s); set(gca,'xscale','log');

%% Partial penetration
Q = well(1).Q; kh = HK(end); kv= VK(end);
b    = abs(well(1).z(2)); d    = abs(well(1).z(1));
bacc = b;                 dacc = d;
D = sum(gr.dz(:));

Ds =0;
for n=1:100
    Ds = Ds + besselk(0,n*pi*rw/D)/n^2 .* (sin(n*pi*b/D)-sin(n*pi*d/D))*(sin(n*pi*bacc/D)-sin(n*pi*dacc/D));
    fprintf('%12.5f\n',Ds);
end
Ds = Ds * Q/(pi*kh)*D/pi^2/((b-d)*(bacc-dacc))*Ds;
%%
figure; hold on;
plot(kDE,'b');
plot(kDS,'r');
plot(kDknown,'k');
legend('kD Elastic','kD Phreatic','kD total');

