% TestFD3DBC  A three dimensional Block centered groundwater flow model
% to 991119, 000409 010517 011003

%Convenient set up of the model, like done with analytic multi-layer models
c =[150;300;500];			% resistance of the aquitards
k =[10;25;35];				% conductivity of the aquitards
z=[0;-10;-40;-50;-80;-100;-150];	% elevation of the layer boundaries (aquitard, aquifer, aquitard, aquifer etc.)
D=abs(diff(z));			% layer thicknesses

K=[D(1)/c(1);k(1);D(3)/c(2);k(2);D(5)/c(3);k(3)];		% layer conductivity also of aquitards

x=logspace(1,log10(2000),20) ; x=[fliplr(-x),x];  % generate x-values of the mesh
y=logspace(1,log10(2000),20)'; y=[flipud(-y);y];  % generate y-values of the mesh
xm=0.5*(x(1:end-1)+x(2:end));
ym=0.5*(y(1:end-1)+y(2:end));
zm=0.5*(z(1:end-1)+z(2:end));

% size of the model mesh
Nx=length(x)-1;
Ny=length(y)-1;
Nz=length(z)-1;

[dum1,dum2,kh]=meshgrid(ones(1,Nx),ones(Ny,1),K); % generate full 3d-grid of conductivity values

FQ=zeros(Ny,Nx,Nz);				% given nodal flows (all zero)
FH=NaN*FQ;							% Fixed head matrix, all NaN's to start with
FH(:,:,1)=0;						% Upper place has fixed head value Phi=0
ii=find(xm==0);
jj=find(ym==0);
kk=find(zm<z(4)&zm>z(5));
FQ(jj,ii,kk)=-2400;
zz=[z(4):-5:z(5)]'

tic;									% remember time
[Phi,Q,QRF,QFF,QLF]=FD3DBC(x,y,z,kh,kh,kh,FH,FQ);			% run het model
toc;									% give time needed for the computation

figure;
ivlak=2;								% contour plane is set to number 2
[cr,dc]=contrange(Phi(:,:,ivlak),25);				% look for a pretty range to contour the heads, using 25 lines
contour(xm,ym,Phi(:,:,ivlak),cr);						% contour the heads
xlabel('x [m]'); ylabel('y [m]');					% add axis labels
title(sprintf('contourlijnen vlak %d  %.1f:%.1f%.1f',ivlak,min(cr),dc,max(cr)));
grid

figure								% another figure
hold on								% prevent destroying for following plot instructions
i=1; patch([x(1),x(end),x(end),x(1)],[z(i),z(i),z(i+1),z(i+1)],[0.7,0.7,0.7]);	% draw first aquitard
i=3; patch([x(1),x(end),x(end),x(1)],[z(i),z(i),z(i+1),z(i+1)],[0.7,0.7,0.7]);	% draw second aquitard
i=5; patch([x(1),x(end),x(end),x(1)],[z(i),z(i),z(i+1),z(i+1)],[0.7,0.7,0.7]);	% draw third aquitard
[cr,dc]=contrange(Phi(20,:,:),50);					% look for a suitable contour set to plot the Psi's
contour(xm,zm,squeeze(Phi(20,:,:))',cr);				% contour them in the vertical plane (y=0)
set(gca,'ylim',[min(z),max(z)]);						% fix the vertical axis scale
xlabel('x [m]'); ylabel('z [mNAP]'); title(sprintf('head contours at %.1f:%.1f:%.1f',cr(1),dc,cr(end)));

disp(['  N.B.: Gebruik zoom xon om in te zoomen op de xas met behoud van de z as']);



% particle tracking
por=0.35*ones(size(kh));
tmax=1e6;
if 0
[XM,YM,ZM]=meshgrid(xm,ym,zm);
plist=[XM(:),YM(:),ZM(:),zeros(size(XM(:)))];
end

if 0
   r=25; theta=[0:15:360]';
   xx=r*cos(theta)/pi; yy=r*sin(theta)/pi; ; tt=zeros(size(theta));
   plist=[];
   for j=1:length(zz)
      plist=[plist;[xx yy zz(j)*ones(size(theta)) tt]];
   end
end

zp=[zm(1):-5:zm(end)]';

plist=[-100*ones(size(zp)) zeros(size(zp)) zp zeros(size(zp));...
   	+100*ones(size(zp)) zeros(size(zp)) zp zeros(size(zp))];



[ptcl]=modpath3(x,y,z,por,QRF,QFF,QLF,plist,tmax)
figure
hold on
for i=1:length(ptcl)
   plot3(ptcl(i).pth(2,:),ptcl(i).pth(3,:),ptcl(i).pth(4,:),'r');
end
