% radmodtest
% thiem
r=logspace(-1,3,20);
y=[0:5:5]';

Nr=length(r);
Ny=length(y);

KR=33.4;
kr=KR*ones(Ny-1,Nr-1);
ky=kr/5;

FH=zeros(Ny,Nr)*NaN;
FQ=zeros(Ny,Nr);
FH(:,end)=-1.35;
FH(:,1)=-5.6;

[Phi,Q]=radmod(r,y,kr,ky,FH,FQ);
PhiRange=[min(Phi(:)):0.1:max(Phi(:))];
contour(r,y,Phi,PhiRange);
Q0=sum(Q(:,1))
Phi0=Phi(1,end)+Q0/(2*pi*KR*(y(end)-y(1)))*log(r(end)/r(1))
semilogx(r,Phi)
%set(gca,'xlim',[0,50],'ylim',[0,50]);