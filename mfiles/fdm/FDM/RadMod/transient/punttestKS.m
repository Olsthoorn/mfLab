% punttest: model run om te zien wat gemeten wordt bij een punttest op een kleine
% peilbuis met open boorgat
% TO 051124

% this version shows effect of different conductivities etc
close all

spd=86400; % seconds per day

%close all

theta=0.67;

t=logspace(-1,3,41); % seconds

rw=0.02; zw=3.0;  % radius and depth of piezometer (flat open bottom)

r=logspace(-3,2,51);                     % r from 1 mm to 10 m
z=fliplr(sort([-3:0.1:5,...
  zw+[-0.15 -0.07 -0.04 -0.02 -0.01 -0.005 0.005 0.01 0.02 0.04 0.07 0.15]]))'; % Grid refined at pipe bottom

Nr=length(r); Nz=length(z); Nt=length(t);

rm=0.5*(r(1:end-1)+r(2:end));
zm=0.5*(z(1:end-1)+z(2:end)); 


lg=[]; clr=['brkgycbrkgycbbrkgycbbrkgycbbrkgycb'];

thisCase=1;

switch thisCase
case 1 % Effect of conductivity
    krCase=[   2    4    7   10   13   15   20];
    kzCase=[   2    4    7   10   13   15   20];
    SsCase=[1e-5 1e-5 1e-5 1e-5 1e-5 1e-5 1e-5];
case 2 % Effect of storage coefficient
    krCase=[  10   10   10   10   10   10   10];
    kzCase=[  10   10   10   10   10   10   10];
    SsCase=[1e-3 1e-4 5e-5 1e-5 5e-6 1e-6 1e-7];
case 3 % anisotropy
    krCase=[  10   10    10   10   10   10   10];
    kzCase=[   2    4     7   10   13   15   20];
    SsCase=[1e-5 1e-5  1e-5 1e-5 1e-5 1e-5 1e-5];
case 4  % Keep overall conductivity constant at 10 m/d
    krCase=[  10   10    10   10   10   10   10];
    kzCase=[   2    4     7   10   13   15   20];
    krCase=sqrt((krCase.^3)./kzCase);
    SsCase=[1e-5 1e-5  1e-5 1e-5 1e-5 1e-5 1e-5];
case 5  % Keep overall conductivity constant at 10 m/d
    krCase=[  10   10    10   10   10   10   10];
    kzCase=[   2    4     7   10   13   15   20];
    krCase=(krCase.^2)./kzCase;
    SsCase=[1e-5 1e-5  1e-5 1e-5 1e-5 1e-5 1e-5];
case 6  % Keep overall conductivity constant at 10 m/d
    krCase=[  10   10    10   10   10   10   10];
    kzCase=[ 0.3     1     3   10   30   100  300];
    krCase=sqrt((krCase.^3)./kzCase);
    SsCase=[1e-5 1e-5  1e-5 1e-5 1e-5 1e-5 1e-5];
case 7  % Keep overall conductivity constant at 10 m/d
    krCase=[  10   10    10   10   10   10   10];
    kzCase=[ 0.3     1     3   10   30   100  300];
    krCase=(krCase.^2)./kzCase;
    SsCase=[1e-5 1e-5  1e-5 1e-5 1e-5 1e-5 1e-5];
case 8
    figure
    Q=1*(pi*rw^2)/200;         % from pipe volume infiltrates in about 200 s during falling head
    k=1e-4; Ss=1e-5;
    rp=rw/4; 
    tau=logspace(-6,2,81),semilogx(tau,Q./((4*pi*k)*rp).*erfc(sqrt(rp.^2*Ss./(4*k*tau))));
    grid on; xlabel('time [s]'); ylable('head change [m]'); 
    title('Analytic solution of sudden continous injection');
end

nCase=length(krCase);
lg=[];

for iCase=1:nCase
    
kr=krCase(iCase)*ones(Nz-1,Nr-1)/spd;
kz=kzCase(iCase)*ones(Nz-1,Nr-1)/spd; % m/s
Ss=SsCase(iCase)*ones(Nz-1,Nr-1);

FH=zeros(Nz,Nr)*NaN; FQ=zeros(Nz,Nr); Phi0=zeros(Nz,Nr);

% find cells inside the pipe
Ir=find(rm<=rw);
Jw=find(zm>=zw);

kr(Jw,Ir)=1e4; kz(Jw,Ir)=1e4;                   % cells in pipe huge conductivity
kr(Jw,Ir(end)+1)=0; kz(Jw,Ir(end)+1)=0;         % cells in pipe wall conductivity = 0

%Storage must be as follows: S=Ss*Dh=1 therfore Ss=1/Dh with Dh=0.1--> Ss=10 applied only in topcell
Ss(Jw(1),Ir)=1/(z(1)-z(2));

J=find(z>=zw); I=find(r<=rw); Phi0(J,I)=1;  % intitial heads inside pipe

[Phi,Q]=radmodt(r,z,t,kr,kz,Ss,FH,FQ,Phi0,theta);

j=find(z>=zw); j=j(end);    % z well bottom

P=permute(Phi,[3,2,1]);
%plot(t,P(:,1,j),clr(iCase)); hold on;      % plot phi(t) points at z=zw
%semilogx(t,P(:,1,j),clr(iCase)); hold on;      % plot phi(t) points at z=zw
semilogy(t,P(:,1,j),clr(iCase)); hold on;      % plot phi(t) points at z=zw
lg=[lg,{sprintf('kr=%.0f,kz=%.0f,Ss=%.1g',krCase(iCase),kzCase(iCase),SsCase(iCase))}];

end % case
legend(lg);
title('head at z=lower end of pipe vs time'); xlabel('time [s]'), ylabel('head [m]');
