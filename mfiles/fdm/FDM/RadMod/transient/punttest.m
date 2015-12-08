% punttest: model run om te zien wat gemeten wordt bij een punttest op een kleine
% peilbuis met open boorgat
% TO 051124

spd=86400; % seconds per day

%close all

theta=0.67;

t=logspace(-5,3,81); % seconds

rw=0.02; zw=3.0;  % radius and depth of piezometer (flat open bottom)

r=logspace(-3,2,51);                     % r from 1 mm to 10 m
z=fliplr(sort([-3:0.1:5,...
  zw+[-0.15 -0.07 -0.04 -0.02 -0.01 -0.005 0.005 0.01 0.02 0.04 0.07 0.15]]))'; % Grid refined at pipe bottom

Nr=length(r); Nz=length(z); Nt=length(t);

rm=0.5*(r(1:end-1)+r(2:end));
zm=0.5*(z(1:end-1)+z(2:end)); 
   
KR=10/spd;      kr=KR*ones(Nz-1,Nr-1); kz=kr; % m/s
SS=1e-5;        Ss=SS*ones(Nz-1,Nr-1);

FH=zeros(Nz,Nr)*NaN; FQ=zeros(Nz,Nr); Phi0=zeros(Nz,Nr);

% find cells inside the pipe
Ir=find(rm<=rw);
Jw=find(zm>=zw);

kr(Jw,Ir)=KR*1e6; kz(Jw,Ir)=KR*1e6;             % cells in pipe huge conductivity
kr(Jw,Ir(end)+1)=0; kz(Jw,Ir(end)+1)=0;         % cells in pipe wall conductivity = 0

%Storage must be as follows: S=Ss*Dh=1 therfore Ss=1/Dh with Dh=0.1--> Ss=10 applied only in topcell
Ss(Jw(1),Ir)=1/(z(1)-z(2));

J=find(z>=zw); I=find(r<=rw); Phi0(J,I)=1;  % intitial heads inside pipe

[Phi,Q]=radmodt(r,z,t,kr,kz,Ss,FH,FQ,Phi0,theta);

j=find(z>=zw); j=j(end);    % z well bottom

P=permute(Phi,[3,2,1]);

lg=[]; clr=['brkgycbrkgycbbrkgycbbrkgycbbrkgycb']; ii=0;
for i=[14 15 18 21 25 28 31]
    ii=ii+1;
    semilogx(t,P(:,i,j),clr(ii)); hold on;      % plot phi(t) points at z=zw
    lg=[lg,{sprintf('r=%.2f',r(i))}];
end
legend(lg);
title('head over time for points at different r at z=lower end of pipe'); xlabel('time [s]'), ylabel('head [m]');

figure
lg=[]; clr=['brkgycbrkgycbbrkgycbbrkgycbbrkgycb']; ii=0;
for j=[27 30 31 33 35 36 38 40 43]
    ii=ii+1;
    semilogx(t,P(:,1,j),clr(ii)); hold on;      % plot phi(t) points at z=zw
    lg=[lg,{sprintf('z=%.2f',z(j))}];
end
legend(lg);
title('head over time for points at different z at r=0'); xlabel('time [s]'), ylabel('head [m]');

% Contour for t=T
figure
T=20; it=find(t<=T); it=it(end);    % locate t=T in time series
PhiRange=round([min(Phi(:)):(max(Phi(:))-min(Phi(:)))/30:max(Phi(:))]*100)/100;
[cs,h]=contour(r,z,Phi(:,:,it(end)),PhiRange);
set(gca,'xlim',[0, 0.5], 'ylim', [2.75 3.25]);
title(sprintf('contours of head at time %5.1f [s]',it(end))); xlabel('r [m]'); ylabel('z [m]');
clabel(cs,h,'manual');

% plot head for different times and z (r=0)
figure
for k=1:length(t)
plot(Phi(:,1,k),z); hold on
end
title('head at r=0 over time'); xlabel('head [m]'), ylabel('z [m]');

% balance of the model
pi*rw^2*max(Ss(:))*(z(1)-z(2))*(Phi(1,1,1)-Phi(1,1,end))
pi*r(end)^2*abs(z(end)-z(1))*min(Ss(:))*(Phi(end,end,end)-Phi(end,end,1))
