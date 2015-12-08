% TestThesis2  A three dimensional Block centered groundwater flow model PhD thesis p107
% to 011004

%Convenient set up of the model, like done with analytic multi-layer models
por=0.35;
x=[0,5,15,25,50,100:100:1000];
y=[0,5,15,25,50,100:100:1500];
z=[0,-5,-40,-50,-100];

D=abs(diff(z));			% layer thicknesses

xm=0.5*(x(1:end-1)+x(2:end));
ym=0.5*(y(1:end-1)+y(2:end));
zm=0.5*(z(1:end-1)+z(2:end));

% size of the model mesh
Nx=length(x)-1;
Ny=length(y)-1;
Nz=length(z)-1;

kh=[0.1, 30, 0.02, 50];  [dum,dum,kh]=meshgrid(ones(1,Nx),ones(1,Ny),kh);
kv=kh; kv(:,:,1)=30;			% vertical conductance of first layer is 30 m/d.

% hole in the aquitarad
I=find(xm>700  & xm<1000); J=find(ym>1200 & ym<1500); kv(J,I,3)=30;

%hole below river
I=find(xm>0  & xm<200); J=find(ym>0 & ym<200); kv(J,I,3)=30;



Por=por*ones(Ny,Nx,Nz);

FQ=  zeros(Ny,Nx,Nz);				% given nodal flows (all zero)
FQ(:,:,1)=0.001*abs(diff(y(:)))*diff(x);		% precipitation on top layer 1 mm/d

FH=  NaN*FQ;							% Fixed head matrix, all NaN's to start with
FH(:,1,1)=0.5.*ym./1000;			% river head drops 0.5 m/km

tic;									% remember time
[Phi,Q,QRF,QFF,QLF]=FD3DBC(x,y,z,kh,kh,kv,FH,FQ);			% run het model
toc;									% give time needed for the computation

figure;
for iz=1:length(zm)
   surf(Phi(:,:,iz)); hold on
end

% particle tracking
times=[3650:3650:200000];
theta=[0:15:360]'/pi;
xp= 700+250*sin(theta);
yp=1100+250*cos(theta);
zp= -5*ones(size(theta));
dp=zeros(size(theta));

plist=[xp yp zp dp];

[ptcl]=modpath3(x,y,z,Por,QRF,QFF,QLF,plist,times)
close all
for i=1:length(ptcl)
   plot3(ptcl(i).pth(2,:),ptcl(i).pth(3,:),ptcl(i).pth(4,:),'r');
   hold on
end

for i=1:length(ptcl)
%   plot3(ptcl(i).tpth(2,:),ptcl(i).tpth(3,:),ptcl(i).tpth(4,:),'r'); hold on
   for j=1:length(ptcl(i).tpth(1,:))
      plot3(ptcl(i).tpth(2,j),ptcl(i).tpth(3,j),ptcl(i).tpth(4,j),'bo');
   end
end
set(gca,'xlim',[0 1000],'ylim',[0 1500],'zlim',[-100 0]);
xlabel('x'); ylabel('y'); zlabel('z');
