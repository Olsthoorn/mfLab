x=[0:100:2500];  Nx=length(x);
y=[0;-25;-50;-75;-80;-100]; Ny=length(y);

k1=0.001; k2=20; k3=0.001; k4=50; k5=50;
kx=[k1,k2,k3,k4,k5]'*ones(1,length(x)-1);
ky=kx;

FQ=zeros(Ny,Nx);
FH=NaN*FQ;
FH(1,:)=0;
FH([2,3,4],1)=[-1,-1,-3]';


[Phi,Q]=flatmeshctrd(x,y,kx,ky,FH,FQ);
hold on

patch([x(1),x(end),x(end),x(1)],[y(2),y(2),y(1),y(1)],[0.85,0.85,0.85]);
patch([x(1),x(end),x(end),x(1)],[y(4),y(4),y(3),y(3)],[0.85,0.85,0.85]);
contour(x,y,Phi,[-3:0.01:0]);


x=logspace(0,3,15);
y=[0,-5,-5-logspace(0,2,10)]';
Nx=length(x);
Ny=length(y);
FQ=zeros(Ny,Nx);
FH=NaN*FQ; FH(1,:)=0; FH([2,3],1)=-5;
kx=20*ones(Ny-1,Nx-1); kx(1,:)=0.1; ky=kx;
[Phi,Q]=flatmeshctrd(x,y,kx,ky,FH,FQ);
close all
patch([x(1),x(end),x(end),x(1)],[y(2),y(2),y(1),y(1)],[0.9,0,9,0.9]);
hold on
contour(x,y,Phi,100);
set(gca,'ylim',[y(end),y(1)]);