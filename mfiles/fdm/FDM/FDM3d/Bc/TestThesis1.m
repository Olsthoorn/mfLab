% TestThesis1  A three dimensional Block centered groundwater flow model PhD thesis p105ff
% to 01104

%Convenient set up of the model, like done with analytic multi-layer models
k =   10;				% conductivity of the aquitards
por=0.35;
x=[0,5,15,25,50,100:100:1000];
y=[0,5,15,25,50,100:100:1000];
z=[0,-5,-10,-25,-50,-100];

D=abs(diff(z));			% layer thicknesses

xm=0.5*(x(1:end-1)+x(2:end));
ym=0.5*(y(1:end-1)+y(2:end));
zm=0.5*(z(1:end-1)+z(2:end));

% size of the model mesh
Nx=length(x)-1;
Ny=length(y)-1;
Nz=length(z)-1;

kh =k  *ones(Ny,Nx,Nz);
Por=por*ones(Ny,Nx,Nz);
FQ=  zeros(Ny,Nx,Nz);				% given nodal flows (all zero)
FH=  NaN*FQ;							% Fixed head matrix, all NaN's to start with
FH(:,end,:)=0;							% yz plane at x=1000, fixed head = 0
FQ(1,1,1)  =-300;						% Extraction in corner of model

tic;									% remember time
[Phi,Q,QRF,QFF,QLF]=FD3DBC(x,y,z,kh,kh,kh,FH,FQ);			% run het model
toc;									% give time needed for the computation

figure;
for iz=1:length(zm)
   surf(Phi(:,:,iz)); hold on
end

% particle tracking
times=[10000:10000:200000];
theta=[0:15:360]'/pi;
xp=900*ones(size(theta));
yp=500+200*cos(theta);
zp=-50+ 30*sin(theta);
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
set(gca,'xlim',[0 1000],'ylim',[0 1000],'zlim',[-100 0]);

