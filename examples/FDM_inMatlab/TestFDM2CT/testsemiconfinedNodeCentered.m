% This mfile runs a test of the finite element model fdm2t implemented
% directly in Matlab. The test is applied on axially symmetric flow.
% First the steady state variant, then the transient version
% Semi-confined flow (Glee and Hantush)
% TO 120110
%

close all
clear

%% grid
y=[-1 0 2];
r=[0 logspace(-2,5,31)]; [r,y,rm,ym,dr,dy,Nr,Ny]=modelsize(r,y);
R=rm(end);

%% Paremetrs
c=1000;
kh=[0;100];
kv=[Inf;Inf];
Q=2400;
S=[0.1;0.001];
lambda=sqrt(kh(2)*dy(2)*c);

%% Numerical boundary conditions
FH=NaN(  Ny,Nr); FH(1,:)=0;
FQ=zeros(Ny,Nr);
FQ(2,1)=Q;

%% Test steady state solution
sSS=Q/(2*pi*kh(2)*dy(2))*besselk(0,rm/lambda);           % Dupuit analytically
phiSSc=fdm2c(r,y,kh,c,kv,FH,FQ,'axial');   % Dupuit numerically
phiSS=fdm2 (r,y,kh,[0.5*dy(1)/c;Inf],FH,FQ,'axial');

figure; hold on; xlabel('r [m]'); ylabel('head [m]'); 
title('Dupuit, Blue:analytic, Red:numeric');

plot(rm,sSS   ,'b');   % analytical
plot(rm,phiSSc ,'r+');  % numerical
plot(rm,phiSS  ,'go');
set(gca,'xscale','log'); grid on;

% Seems to be perfectly correct, however

set(gca,'yscale','log');

%%

% the numerical solution deviates for very small values of the drawdown.
% This is probably due to the solver. If it breaks off at some small value,
% we will have an error for this and dlower values. This cannot be sovled
% unless we solve in an exact fashion.

fh=NaN(Ny+1,Nr+1);
fq=zeros(Ny+1,Nr+1); fq(:,1)=Q/2;

[Phi,Q,Psi,Qx,Qy]=fdm2NodeCentered(r,y,kh,kv,fh,fq,2);

return
%% test transient solution

t=logspace(-3,2,51)'; % time

s=Q/(4*pi*T)*expint(S./(4*T*t)*(rm.^2));  % Theis analytical solution

% numerical boundary conditions
IH=s(1,:);          % initial head
FH=NaN(size(FH));   % fixed heads

Phi = fdm2t(r,y,t,T,T,S,IH,FH,FQ,'axial'); Phi=permute(Phi,[3,2,1]); % Theis numerically

figure; hold on; xlabel('time [d]'); ylabel('head [m]'); title('Theis, Blue:analytic, Red:numeric');
set(gca,'xscale','log','xgrid','on','ygrid','on');

plot(t,s,   'b'); % Theis analytically as function of time
plot(t,Phi,'r+'); % Theis numerically  as function of time

figure; hold on; xlabel('distance [m]'); ylabel('head [m]'); title('Theis, Blue:analytic, Red:numeric');
set(gca,'xscale','log','xgrid','on','ygrid','on');

plot(rm,s,   'b'); % Theis analytically as function of r
plot(rm,Phi,'r+'); % Theis numerically  as function of r

%% Everything seems to look great, however ...

set(gca,'yscale','log');

% shows that values near zero tend to be inaccurate.
% This can be solved by taking many more points along the r-axis.
% the model seems sensitive to this for small drawdowns.
% take the number of points beyond 2000 and the entire line
% will be accurately computed by the numeric model.

%% Test semiconfined 2



