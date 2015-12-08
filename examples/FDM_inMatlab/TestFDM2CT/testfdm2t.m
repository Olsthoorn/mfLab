% This mfile runs a test of the finite element model fdm2t implemented
% directly in Matlab. The test is applied on axially symmetric flow.
% First the steady state variant, then the transient version
% TO 120110
%
close all
clear

%% Paremetrs
T=100;
Q=2400;
S=0.1;

%% grid
y=[-1 0];
r=logspace(-2,5,71); [r,y,rm,ym,dr,dy,Nr,Ny]=modelsize(r,y);
R=rm(end);

%% Numerical boundary conditions
FH=NaN(size(rm)); FH(end)=0;
FQ=zeros(size(rm));
FQ(1)=Q;


%% Test steady state solution

sSS=Q/(2*pi*T)*log(R./rm);           % Dupuit analytically
phiSS=fdm2(r,y,T,T,FH,FQ,'axial');   % Dupuit numerically

figure; hold on; xlabel('r [m]'); ylabel('head [m]'); 
title('Dupuit, Blue:analytic, Red:numeric');

plot(rm,sSS   ,'b');   % analytical
plot(rm,phiSS ,'r+');  % numerical

set(gca,'xscale','log'); grid on;
%% test transient solution

t=logspace(-3,2,51)'; % time

s=Q/(4*pi*T)*expint(S./(4*T*t)*(rm.^2));  % Theis analytical solution

% numerical boundary conditions
IH=s(1,:);          % initial head
FH=NaN(size(FH));   % fixed heads

Phi = fdm2t(r,y,t,T,T,S,IH,FH,FQ,'axial'); Phi=permute(Phi,[3,2,1]); % Theis numerically

figure; hold on; xlabel('time [d]'); ylabel('head [m]'); title('Theis, Blue:analytic, Red:numeric');
set(gca,'xscale','log','xgrid','on','ygrid','on');

plot(t,s/(Q/(4*pi*T)),'b'); % Theis analytically as function of time
plot(t,Phi/(Q/4*pi*T),'r+'); % Theis numerically  as function of time

figure; hold on; xlabel('distance [m]'); ylabel('head [m]'); title('Theis, Blue:analytic, Red:numeric');
set(gca,'xscale','log','xgrid','on','ygrid','on');

plot(rm,s/(Q/(4*pi*T)),   'b'); % Theis analytically as function of r
plot(rm,Phi/(Q/(4*pi*T)),'r+'); % Theis numerically  as function of r

%% Everything seems to look great, however ...

%set(gca,'yscale','log');

% shows that values near zero tend to be inaccurate. Although this does not
% matter in practice, it is a fundamental problem, that I have not been
% able to solve. It may be a numerical problem with the solver, showing too
% little accuracy for low end values of the drawdown.

