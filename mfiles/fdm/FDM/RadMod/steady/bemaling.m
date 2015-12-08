%Bemaling distributie Ed De Greef 050202

R=10.6 % equivalent radisu of pit
L=12;  % length of sheet piling
z1=-4.5% base of Holocene
h1=-2; % initial head
h2=-4; % final head
kh=40; % conductivity
kv=40; % conductivity
c =500;% resistance
dy=0.5;% thickness of model layers
kc=0.5*dy/c % we put FH in lowest row of holocene

r=logspace(-1,4,60);   rm=0.5*(r(1:end-1)+r(2:end));
y=[0:-dy:-40]';       ym=0.5*(y(1:end-1)+y(2:end));


idw=find(abs(rm-R)==min(abs(rm-R)));
jdw=find(ym>-12);
jtop=find(ym>-4.5);
jaqu=find(ym<-4.5);
iwell=idw-1;
jwell=find(ym<-4.5 & ym>-12);

k=zeros(length(ym),length(rm));
k(jaqu,:)=40;
k(jtop,:)=kc;
%k(jdw,idw)=1e-4;
FH=zeros(size(k))*NaN;
FQ=zeros(size(k));
FH(jtop,:)=-2;
FH(jwell,iwell)=-4;

[Phi,Q,QRF,QFF]=radblockctrd(r,y,k,k,FH,FQ);
[Psi]=streamf(r,y,QRF,QFF);

close all
contour(rm,ym,Phi,30);
hold on
contour(r,y,Psi,30)
set(gca,'xscale','log');



Qextr=sum(Q(jwell,iwell))/24;

