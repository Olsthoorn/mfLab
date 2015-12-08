% for the modeflow keynote in 2001
% comparing meshctrd with blockcntrd and analytical
% as in the thesis
% TO 001019

close all
n=0.001;
L=600;

x=logspace(0,log10(L),10);
x=[0:L/30:L];
y=[0;-1;-2;-3];
% Block centerd model
x1=[-0.01,0,x,L+0.01];										% add small cells to ends to fix boundaries

xm=0.5*(x1(1:end-1)+x1(2:end));							% cell centers
y1=[100;0];
dx1=diff(x1);													% cell width
dy1=abs(diff(y1));
Nx=length(dx1);												% number of cells in a row
Ny=length(dy1);												% number of cells in a column
k1=30*ones(Ny,Nx);											% cell conductivities
FH1=NaN*ones(Ny,Nx); FH1(:,1)=0; FH1(:,end)=0;		% fixed heads
FQ1=(dy1*dx1)*n;												% fixed flows (precipitation flow per cell)
[Phi1,Q1]=flatblockctrd(x1,y1,k1,k1,FH1,FQ1);		% solution for cells


% Mesh centerd grid
x2=x;
y2=y;
dx2=diff(x2);
dy2=abs(diff(y2));
Nx2=length(x2);
Ny2=length(y2);
k2=30*ones(Ny2-1,Nx2-1);									% element conductivities
FH2=NaN*ones(Ny2,Nx2); FH2(:,1)=0; FH2(:,end)=0;	% fixed heads
FQ2=n*([dy2*dx2,0]+[0,dy2*dx2])/4;
FQ2=[FQ2;FQ2];	% fixed flows,
[Phi2,Q2]=flatmeshctrd(x2,y2,k2,k2,FH2,FQ2);			% solution

%analytisch
ka=30
xa=[0:L/50:L];
h=-n*((xa-L/2).^2-(L/2).^2)/(2*ka);


plot(xm,Phi1,'r'); hold on
plot(x2,Phi2,'b');
plot(xa,h   ,'g');