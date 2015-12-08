% radmodt test
% simulating delayed yield and generating Boulton/Neuman type curves
% TO 990522 051221

close all


r=0.5*logspace(-1,3,41);                rm=0.5*(r(1:end-1)+r(2:end));
z=[0:-0.2:-8]'; dz=abs(diff(z));    zm=0.5*(z(1:end-1)+z(2:end));
Nr=length(r);
Nz=length(z);

ZSB=-8; ZST=1; % screen
Q=200;

J=find(z>=ZSB & z<=ZST); q=Q/(length(J)-1);

Sy=0.25; SS=0.00001;
kr=10  *ones(Nz-1,Nr-1); ky=kr;
ss=SS*ones(Nz-1,Nr-1); ss(1,:)=Sy./dz(1,:);

t=logspace(-6,0,61)';

Phi0=zeros(Nz,Nr);											% initialize with theis(r,t(1));
FQ=zeros(Nz,Nr); FQ(J,1)=q; FQ([J(1),J(end)],1)=q/2;
FH  =zeros(Nz,Nr)*NaN;										% no fixed heads

theta=1;
[Phi,Q]=radmodt(r,z,t,kr,ky,ss,FH,FQ,Phi0,theta);
M=max(Phi(:)); m=min(Phi(:)); N=50; C=[m:(M-m)/N:M];
clear H;
for i=1:size(Phi,3)-1
    contour(r,z,Phi(:,:,i+1),C);
    H(i)=getframe(gcf);
end
movie(H);
figure
semilogx(t,shiftdim(Phi(end,1,:),1));
figure
loglog(t,shiftdim(Phi(end,1,:),1));
