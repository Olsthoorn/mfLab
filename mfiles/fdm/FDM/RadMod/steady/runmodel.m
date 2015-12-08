% Bouwput waterleidingkanaal
% TO 990522

% Bouwput, cirkelvormig, R=8.5 m
% kD=1600; D=44;
% Damwand steekt 6 m in het wvp
% c top = 100d
% we gaan uit van putten met 5 m filter aan binnenzijde van de damwand en rond tip damwand

close all

kD=1600
D =44
c =100
hTop=-1.35
hPutten=-6.35;
KR=kD/D
KY=KR/2;
lambda=sqrt(kD*c)
Rbp=8.5;
yTop=0;
yBot=yTop-D;
yDw=-6
yPut=-5
PhiBp=-5.6;

% fijn grid r
r1=Rbp+0.05-fliplr(logspace(-1,log10(Rbp-0.05),15));
r2=Rbp-0.05+logspace(-1,log10(1000-Rbp+0.05),30);
r=[r1,r2];

% fijn grid y
y0=1;
y1=fliplr([0,logspace(-1,log10(-yDw),11)]+yDw);
y2=yDw-logspace(-1,log10(-yBot+yDw),25);
y=[y0,y1,y2]';

Nr=length(r);
Ny=length(y);

rc=(r(1:end-1)+r(2:end))/2;
yc=(y(1:end-1)+y(2:end))/2;
i=(rc<Rbp+0.04 & rc>Rbp-0.04);
j=(yc>yDw-0.01& yc<=-0.01);
I=find(ones(size(yc))*i & j*ones(size(rc)));


kr=KR*ones(Ny-1,Nr-1);
ky=KY*ones(Ny-1,Nr-1);
kr(1,:)=(y(1)-y(2))/c;
ky(1,:)=(y(1)-y(2))/c;
kr(I)=0;
ky(I)=0;
ky(1,find(rc<Rbp+0.05))=0;

i=(r<Rbp-0.01 & r>Rbp-0.06);
j=(y>=yPut & y<=0);
I=find(ones(size(y))*i & j*ones(size(r)));

FH=zeros(Ny,Nr)*NaN; FH(1,:)=hTop; FH(I)=hPutten;
FQ=zeros(Ny,Nr);

[Phi,Q]=radmod(r,y,kr,ky,FH,FQ);
Qo=sum(sum(Q(2:end,:)))
Qt=sum(Q(1,:))

s=sprintf('stijghoogte radiaal model, en De Glee, kr=%3.1fm2/d kz=%3.1fm2/d Q=%4.0fm3/d hTop=%2.2fm hPutten=%2.2fm\n',...
           KR,KY,Qo,hTop,hPutten);

PhiRange=[min(Phi(:)):0.1:max(Phi(:))];
figure; contour(r,y,Phi,PhiRange);
axis equal;
set(gca,'xlim',[0,100],'ylim',[-D 0]); grid on
title(s); xlabel('r [m]'); ylabel('z [m]');


Phi0=Qo/(2*pi*kD)*besselk(0,r/lambda)+hTop;
figure; plot(r,Phi,r,Phi0,'k'); grid on
title(s);
xlabel('r [m]'); ylabel('z [m]');

% Semipervious wall resistance
% De Glee
PhiBpRand=Qo/(2*pi*kD)*besselk(0,Rbp/lambda);
v=Qo/(2*pi*Rbp*D);
Dh=PhiBp-PhiBpRand;
c=Dh/v



