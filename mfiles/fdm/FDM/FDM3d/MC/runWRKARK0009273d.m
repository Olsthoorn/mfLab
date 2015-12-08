% testmc3d  A three dimensional mesh-centerd groundwater flow model
% to 991119, 000409
clear
close all

%Convenient set up of the model, like done with analytic multi-layer models
c =[500];			% resistance of the aquitards
k =[55;];				% conductivity of the aquitards
z=[0;-3.5;-55];	% elevation of the layer boundaries (aquitard, aquifer, aquitard, aquifer etc.)
D=abs(diff(z));			% layer thicknesses

K=[D(1)/c(1);k(1)];		% layer conductivity also of aquitards
% ====Building pit size ====
B =20;				% Length sleuf in x-direction
L =40;				% Length of buidling pit in y-direction
h =-5;				% Head at bottom of building pit to be maintained

y1=[10,15,17.5,18.5,19,19.5,19.75,20,20.25,20.5,20.75,21,22,24,27,30,40,50,75,logspace(2,log10(4000),10)];

x1=[0,2,5,10,15,20,22.5,25,27,28,29,29.5,30,30.5,31,32,32.5,33,33.5,34,35,37.5,40,45,50,60,75,...
      logspace(2,log10(4000),10)];
x2=[0.5,1,2,4,7.5,10,15,25,50,75,100,125,130,135,137,138,139,139.5,140,141,143,145,150,160,175,...
      logspace(log10(200),log10(4000),10)];

y=flipud(y1');  % generate x-values of the mesh
x=[fliplr(-x2),x1];  % generate y-values of the mesh

% size of the model mesh
Nx=length(x);
Ny=length(y);
Nz=length(z);

[dum1,dum2,kh]=meshgrid(ones(1,Nx-1),ones(Ny-1,1),K); % generate full 3d-grid of conductivity values
clear dum1 dum2




xm=0.5*(x(1:end-1)+x(2:end));
ym=0.5*(y(1:end-1)+y(2:end));
I1=find(xm<33 & xm>  30); kh(:,I1,1)=3.5/0.001;
I2=find(xm<0  & xm>-140); kh(:,I2,1)=3.5/40;

FQ=zeros(Ny,Nx,Nz);				% given nodal flows (all zero)
FH=NaN*FQ;							% Fixed head matrix, all NaN's to start with
FH(:,:,1)=-0.4;						% Upper place has fixed head value Phi=0
I=find(x>=30 & x<= 33);	% Look for the x-nodes within the building pit
J=find(y<=20);
for j=J
   FH(j,I,1)=-2.5;
end

I=find(x>=-140 & x<=0);	  FH(:,I,1)=-0.0;	% Same for y

tic;									% remember time
[Phi,Q]=mc3dfdm(x,y,z,kh,kh,kh,FH,FQ);			% run het model
toc;									% give time needed for the computation

figure;
ivlak=2;								% contour plane is set to number 2
[cr,dc]=contrange(Phi(:,:,ivlak),25);				% look for a pretty range to contour the heads, using 25 lines
contour(x,y,Phi(:,:,ivlak),[-2.5:0.025:0]);						% contour the heads
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
