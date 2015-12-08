% radmodtCC test
% simulating delayed yield and generating Boulton/Neuman type curves
% TO 990522 051221

close all

r=0.4*logspace(-1,2,31);            rm=0.5*(r(1:end-1)+r(2:end)); rm(1,end)=r(1,end);
z=[0:-0.2:-8]'; dz=abs(diff(z));    zm=0.5*(z(1:end-1)+z(2:end)); zm(1,end)=z(1,end);
Nr=length(r);
Nz=length(z);

ZSB=-8; ZST=0; % screen
Q=200;

J=find(zm>=ZSB & zm<=ZST); q=Q/length(J);

Sy=0.25; SS=0.00001;
kr=10*ones(Nz-1,Nr-1); kz=kr;
ss=SS*ones(Nz-1,Nr-1); ss(1,:)=Sy./dz(1,:);
por=0.35*ones(Nz-1,Nr-1);
t=logspace(-6,0,61)';

Phi0=zeros(Nz-1,Nr-1);											% initialize with theis(r,t(1));
FQ=zeros(Nz-1,Nr-1); FQ(J,1)=q;
FH  =zeros(Nz-1,Nr-1)*NaN;										% no fixed heads

theta=1;
[Phi,Q,QR,QZ,St]=radmodtCC(r,z,t,kr,kz,ss,FH,FQ,Phi0,theta);
[vr,vz]=radmodtCCQuiver(r,z,t,kr,kz,por,QR,QZ);

M=max(Phi(:)); m=min(Phi(:)); N=50; C=[m:(M-m)/N:M];
clear H;
for it=1:length(diff(t))
    contour(rm,zm,Phi(:,:,it+1),C); hold on
%    quiver(rm,zm,100*vr(:,:,it),100*vz(:,:,it));
    H(i)=getframe(gcf);
end
movie(H);
figure
semilogx(t,shiftdim(Phi(end,1,:),1));
figure
loglog(t,shiftdim(Phi(4,1,:),1),'r'); hold on
loglog(t,shiftdim(Phi(4,13,:),1),'b'); hold on
loglog(t,shiftdim(Phi(21,1,:),1),'r:'); hold on
loglog(t,shiftdim(Phi(21,13,:),1),'b:'); hold on
loglog(t,shiftdim(Phi(end,1,:),1),'r--'); hold on
loglog(t,shiftdim(Phi(end,13,:),1),'b--'); hold on

