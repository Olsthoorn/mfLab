% test_FD3DMC  A three dimensional mesh-centerd groundwater flow model
% to 991119, 000409 TO 001520

%Convenient set up of the model, like done with analytic multi-layer models
c =[150;300;500];			% resistance of the aquitards
k =[10;25;35];				% conductivity of the aquitards
z=[0;-10;-40;-50;-80;-100;-150];	% elevation of the layer boundaries (aquitard, aquifer, aquitard, aquifer etc.)
D=abs(diff(z));			% layer thicknesses

K=[D(1)/c(1);k(1);D(3)/c(2);k(2);D(5)/c(3);k(3)];		% layer conductivity also of aquitards
% ====Building pit size ====
B =20;				% Length of building pit in x-direction
L =40;				% Length of buidling pit in y-direction
h =-5;				% Head at bottom of building pit to be maintained

x=logspace(1,log10(2000),20) ; x=[fliplr(-x),x];  % generate x-values of the mesh
y=logspace(1,log10(2000),20)'; y=[flipud(-y);y];  % generate y-values of the mesh

% size of the model mesh
Nx=length(x);
Ny=length(y);
Nz=length(z);

[dum1,dum2,kh]=meshgrid(ones(1,Nx-1),ones(Ny-1,1),K); % generate full 3d-grid of conductivity values

FQ=zeros(Ny,Nx,Nz);				% given nodal flows (all zero)
FH=NaN*FQ;							% Fixed head matrix, all NaN's to start with
FH(:,:,1)=0;						% Upper place has fixed head value Phi=0
I=find(x>=-B/2 & x<=B/2);		% Look for the x-nodes within the building pit
J=find(y>=-L/2 & y<=L/2);		% Same for y
FH(J,I,2)=h;						% These get the building pit maintained head (in model node layer 2)

tic;									% remember time
[Phi,Q]=FD3DMC(x,y,z,kh,kh,kh,FH,FQ);			% run het model
toc;									% give time needed for the computation

figure;
ivlak=2;								% contour plane is set to number 2
[cr,dc]=contrange(Phi(:,:,ivlak),25);				% look for a pretty range to contour the heads, using 25 lines
contour(x,y,Phi(:,:,ivlak),cr);						% contour the heads
xlabel('x [m]'); ylabel('y [m]');					% add axis labels
title(sprintf('contourlijnen vlak %d  %.1f:%.1f%.1f',ivlak,min(cr),dc,max(cr)));
grid

figure								% another figure
hold on								% prevent destroying for following plot instructions
i=1; patch([x(1),x(end),x(end),x(1)],[z(i),z(i),z(i+1),z(i+1)],[0.7,0.7,0.7]);	% draw first aquitard
i=3; patch([x(1),x(end),x(end),x(1)],[z(i),z(i),z(i+1),z(i+1)],[0.7,0.7,0.7]);	% draw second aquitard
i=5; patch([x(1),x(end),x(end),x(1)],[z(i),z(i),z(i+1),z(i+1)],[0.7,0.7,0.7]);	% draw third aquitard
[cr,dc]=contrange(Phi(20,:,:),50);					% look for a suitable contour set to plot the Psi's
contour(x,z,squeeze(Phi(20,:,:))',cr);				% contour them in the vertical plane (y=0)
set(gca,'ylim',[min(z),max(z)]);						% fix the vertical axis scale
xlabel('x [m]'); ylabel('z [mNAP]'); title(sprintf('head contours at %.1f:%.1f:%.1f',cr(1),dc,cr(end)));

disp(['  N.B.: Gebruik zoom xon om in te zoomen op de xas met behoud van de z as']);
