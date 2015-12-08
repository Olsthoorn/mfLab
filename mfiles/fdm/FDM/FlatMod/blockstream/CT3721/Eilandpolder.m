% test model blkstrm, block centred finite difference model with stream function
% TO 000530, 001026, stroomfunctie ingebouwd.
% TO 050321  Doorsnede model stadsinrichting Student Kooijman voor vak

B=100; % afstand tussen de sloten
mv=0.5; % maaiveld
sp=0.0; % slootpeil
sb=-0.5; % slootbodem
bs=-2.0; % basis
infiltratie=0.001;
Y=[mv,sp,sb,bs];  % laaggrenzen (vallen samen met slootpeil en slootbodem (toeval)

x=[0 0.05 0.5 0.8 0.9 0.95 1 1.1 1.3 1.5:0.5:B/2]; xm=0.5*(x(1:end-1)+x(2:end)); Nx=length(x)-1;						% dx(1)=0 to put boundary at boundary of model
dx=diff(x);
y=[Y(1):-0.05:Y(end)]'; ym=0.5*(y(1:end-1)+y(2:end));  Ny=length(y)-1;

isl=find(xm<1);
jsl=find(ym>-0.5);

i1=find(ym<Y(1) & ym>Y(2));
i2=find(ym<Y(2) & ym>Y(3));
i3=find(ym<Y(3));

XP=[x(1) x(end) x(end) x(1)];
Color=[1,1,0; 1,0.8,0.2; 1,0.6,0.4]; 

k1=2  ;  % 10
k2=0.1; %  3
k3=2  ; %  5

Kx=ones(Ny,Nx);
Kx(i1,:)=k1;
Kx(i2,:)=k2;
Kx(i3,:)=k3;

Ky=Kx/2;  %5


FH=NaN*ones(size(Kx));
FH(jsl,isl)=slootpeil;

FQ=zeros(size(Kx)); FQ(1,:)=infiltratie.*dx;

[Phi,Q,Psi]=blckstrm(x,y,Kx,Ky,FH,FQ);

%close all
figure
for i=1:3; patch(XP,[Y(i),Y(i),Y(i+1),Y(i+1)],Color(i,:)); hold on; end
contour(xm,ym,Phi,25)
hold on
contour(x(2:end-1),y,Psi(:,2:end-1),20);
plot(xm,Phi(1,:))
patch([0,1,1,0],[sb,sb,sp,sp],'b');
patch([0,1,1,0],[sp,sp,mv,mv],'w');
axis equal; axis tight;

title('alle cellen doen mee');
Z=y*ones(size(dx));

f=(Phi-Z(2:end,:))./(Z(1:end-1,:)-Z(2:end,:)); I=find(f<=0); f=max(0.001,min(1,f)); 

[Phi,Q,Psi]=blckstrm(x,y,f.*Kx,Ky,FH,FQ);

h=Phi(1,:); Phi(I)=NaN;
%close all
figure;
for i=1:3; patch(XP,[Y(i),Y(i),Y(i+1),Y(i+1)],Color(i,:)); hold on; end
contour(xm,ym,Phi,25); hold on
contour(x(2:end-1),y,Psi(:,2:end-1),20);
plot(xm,h)
patch([0,1,1,0],[sb,sb,sp,sp],'b');
patch([0,1,1,0],[sp,sp,mv,mv],'w');
%axis equal; axis tight;

title('droge cellen uitgeschakeld, halfdroge cellen in kx gereduceerd');

figure; % alleen opbolling
plot(xm,h)
title('opbolling');
