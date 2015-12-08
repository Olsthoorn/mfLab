% test model blkstrm, block centred finite difference model with stream function
% TO 040905, 001026, stroomfunctie ingebouwd.

close all
n=0.02;
L=8000;
Z=-1000;

x=[0,0:25:L]; Nx=length(x)-1;						% dx(1)=0 to put boundary at boundary of model
y=[0,0,-10:-10:Z,Z]'; Ny=length(y)-1;			% some dy=0 to get proper contours
xm=0.5*(x(1:end-1)+x(2:end));
ym=0.5*(y(1:end-1)+y(2:end));

kx=30*ones(Ny,Nx);
ky=kx;

FH=zeros(Ny,Nx)*NaN; 
FQ=zeros(Ny,Nx);
if 0
FH(1,:)=0;
for i=1:5
    m=rand*10;
    FH(1,:)=FH(1,:)+sin(2*pi*m*xm/L);
end
FQ=zeros(Ny,Nx);
else
    FQ(1,:)=n*diff(x);
    Nx=length(xm);
    R=round(Nx*rand(1,Nx/10)); R=max(1,R); R=min(Nx,R);
    FH(1,R)=rand(1,length(R));
end


[Phi,Q,Psi]=blckstrm(x,y,kx,ky,FH,FQ);
subplot(2,1,1); plot(xm,Phi(1,:));

subplot(2,1,2);
contour(xm,ym,Phi,50)
hold on
contour(x(2:end-1),y,Psi(:,2:end-1),40);

axis equal
axis tight
