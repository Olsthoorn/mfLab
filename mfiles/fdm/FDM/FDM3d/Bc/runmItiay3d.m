% testmc3d  A three dimensional mesh-centerd groundwater flow model
% to 991119, 000409
% to 001128 putten Itiay el Barud

%Convenient set up of the model, like done with analytic multi-layer models

x=[0,logspace(1,log10(200),20),201];						% oplopende celgrootte, extra vlak op x=200;
y=[0,logspace(1,log10(200),20),201]';						% oplopende celgrootte, extra vlak op y=200;
z=[0,-0.1,-14.9,[-15:-5:-300]]';								% Dun bovenvlak, kleilaag en dan zand
xm=0.5*(x(1:end-1)+x(2:end));
ym=0.5*(y(1:end-1)+y(2:end));
zm=0.5*(z(1:end-1)+z(2:end));

% size of the model in cellen
Nx=length(xm);
Ny=length(ym);
Nz=length(zm);

k=30*ones(1,1,Nz); k(1:3)=0.015;
[dum1,dum2,K]=meshgrid(ones(1,Nx),ones(Ny,1),k); % generate full 3d-grid of conductivity values

FQ=zeros(Ny,Nx,Nz);			% given nodal flows (all zero)
FH=NaN*FQ;							% Fixed head matrix, all NaN's to start with
FH(  :,  :,1)=0;					% zm=zm(1) has fixed head value Phi=0
FH(end,  :,:)=0;					% ym=ym(end) idem
FH(  :,end,:)=0;					% xm=xm(end) idem
FH(1,1,[6:14])=-2;				% top screen
FH(1,1,[16:24])=2;				% bottom screen

tic;									% remember time
[Phi,Q,QRF,QFF,QLF]=fd3dbc(x,y,z,K,K,K,FH,FQ);			% run het model
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
i=3; patch([x(1),x(end),x(end),x(1)],[z(i),z(i),z(1),z(1)],[0.7,0.7,0.7]);	% draw first aquitard
[cr,dc]=contrange(Phi(20,:,:),50);					% look for a suitable contour set to plot the Psi's
contour(xm,zm,squeeze(Phi(end,:,:))',cr);				% contour them in the vertical plane (y=0)
set(gca,'ylim',[min(z),max(z)]);						% fix the vertical axis scale
xlabel('x [m]'); ylabel('z [mNAP]'); title(sprintf('head contours at %.1f:%.1f:%.1f',cr(1),dc,cr(end)));
disp(['  N.B.: Gebruik zoom xon om in te zoomen op de xas met behoud van de z as']);

[XM,YM,ZM]=meshgrid(xm,ym,zm);
figure
contourslice(XM,YM,ZM,Phi, [5],[5],[-47.5],50);
view(-20,20); axis vis3d tight                             
daspect([1 1 1]);
box on                                   
camproj perspective; camva(6)                              
colormap(jet)

%Streamlines see matlab demo                  
figure
zp=zm([6:14,16:24]); d=zeros(size(zp));
plist=[xm(1)*ones(size(zp)) y(2) *ones(size(zp)) zp d ; ...
		 x(2) *ones(size(zp)) ym(1)*ones(size(zp)) zp d];
 por=0.35*ones(Ny,Nx,Nz);
 
figure
times=[1:100];
[ptcl]=modpath3(x,y,z,por,QRF,QFF,QLF,plist,times)
for i=1:length(ptcl)
   plot3(ptcl(i).pth(2,:),ptcl(i).pth(3,:),ptcl(i).pth(4,:),'r'); hold on
   for j=1:length(ptcl(i).tpth(1,:))
      plot3(ptcl(i).tpth(2,j),ptcl(i).tpth(3,j),ptcl(i).tpth(4,j),'bo');
   end
end
set(gca,'xlim',[x(1) x(end)],'ylim',[y(1) y(end)],'zlim',[z(end) z(1)]);