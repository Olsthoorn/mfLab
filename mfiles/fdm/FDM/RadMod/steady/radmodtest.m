% radmodtest
% thiem
r=logspace(-1,3,20);
y=[0:1:5]';

Nr=length(r);
Ny=length(y);

KR=33.4;
kr=KR*ones(Ny-1,Nr-1);
ky=kr*100;
ky(1,:)=500/5;
FH=zeros(Ny,Nr)*NaN;
FH(1,:)=0;

FQ=zeros(Ny,Nr);
FQ(2:3,1)=-50;


[Phi,Q]=radmod(r,y,kr,ky,FH,FQ);

close all
PhiRange=[min(Phi(:)):0.1:max(Phi(:))];
contour(r,y,Phi,PhiRange);
Q0=sum(Q(:,1))
Phi0=Phi(1,end)+Q0/(2*pi*KR*(y(end)-y(1)))*log(r(end)/r(1))
semilogx(r,Phi)
legend('1','2','3','4','5');
%set(gca,'xlim',[0,50],'ylim',[0,50]);