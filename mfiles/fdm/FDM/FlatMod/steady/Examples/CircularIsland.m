% Compare anlytical solution for steady-state head with model flatmod
% Problem suggested by Mark Bakker, 060112, because Modflow needs many cells for accurate solution
% Use fixed transmissivity
% TO 060114

% conclusion: Voor R=1000 m 20m cells are necessary for accurate computation,
% plus a small shift in R by 0.4 cell size to get the boundary on its place (prevent shift of entire solution)
% this last thing was done experimentall.
% A much better solution will be obtained by a radially symmetric model, that's for sure, but not always possible.

clear; close all
path(path,'C:\Documents and Settings\POHY01\Desktop\MODELLEN\FDM\flatmod\steady');

L=1100; % Half Size of model
Ro=1000; % Radius of island
D=20;  % cell size
T=100;    % transmissivity
N=0.001;

%Analytic solution

%Numeric solution
x=linspace(-L,L,2*L/D+1); y=fliplr(x)';
xm=0.5*(x(1:end-1)+x(2:end)); ym=0.5*(y(1:end-1)+y(2:end)); [Xm,Ym]=meshgrid(xm,ym);Rm=sqrt(Xm.^2+Ym.^2);
Nx=length(xm); Ny=length(ym);
k=T*ones(Ny,Nx);
FH=NaN*ones(Ny,Nx); FH(find(Rm>Ro-0.4*D))=0;
FQ=N*abs(diff(y))*diff(x);

[Phi,Q,QRF,QFF]=flatblockctrd(x,y,k,k,FH,FQ);
PhiA=N*(Ro^2-Rm.^2)/(4*T); PhiA(find(Rm>Ro))=0% analytic

close all
c=[0:0.1:3];
Nh=round(length(ym)/2);
figure; hold on
contour(xm,ym(1:Nh),Phi(1:Nh,:),c);      % top half of figure numeric model result
contour(xm,ym(Nh:end),PhiA(Nh:end,:),c); % lower half of figure analytic model result


