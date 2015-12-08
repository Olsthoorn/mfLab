% test model blkstrm, block centred finite difference model with stream function
% TO 000530, 001026, stroomfunctie ingebouwd.

x=[0,0:25:1000]; Nx=length(x)-1;						% dx(1)=0 to put boundary at boundary of model
y=[0,0,-10,-10,-100,-100]'; Ny=length(y)-1;			% some dy=0 to get proper contours

kx=30*ones(Ny,Nx); kx(2,:)=0.05;
ky=kx;

FH=zeros(Ny,Nx)*NaN; FH(1,:)=0; FH(4,20)=-1; FH(4,1)=1;
FQ=zeros(Ny,Nx);

[Phi,Q,Psi]=blckstrm(x,y,kx,ky,FH,FQ);

xm=0.5*(x(1:end-1)+x(2:end));
ym=0.5*(y(1:end-1)+y(2:end));

close all
contour(xm,ym,Phi,50)
hold on
contour(x(2:end-1),y,Psi(:,2:end-1),40);
figure;												% contouren stijghoogten en stroomfunctie met verticale cuts

