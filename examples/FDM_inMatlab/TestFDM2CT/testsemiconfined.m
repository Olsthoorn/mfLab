% This mfile runs a test of the finite element model fdm2t implemented
% directly in Matlab. The test is applied on axially symmetric flow.
% First the steady state variant, then the transient version
% Semi-confined flow (Glee and Hantush)
% TO 120110


close all
clear

%% Paremetrs
T=1000;
c=1000;
b=1;
kh=[0;T/b];
kv=[0.5*b/c;Inf];
Q=2400;
lambda=sqrt(T*c);

r1=1e-6*lambda;
r2=1e2 *lambda;
n=100;

r=logspace(log10(r1),log10(r2),n);
y=[-b 0 1];
[r,y,rm,ym,dr,dy,Nr,Ny]=modelsize(r,y);

B=dy;

%% Numerical boundary conditions
IBOUND=ones(Ny  ,Nr  ); IBOUND(1,:)=-1;
ibound=ones(Ny+1,Nr+1); ibound(1,:)=-1;

IH=zeros(  Ny,Nr);           ih=zeros  (Ny+1,Nr+1); fh(1,:)=0;
FQ=zeros(Ny,Nr);             fq=zeros(Ny+1,Nr+1);
FQ(end,1)=Q;                 fq(2:end,1)=Q/2;

%N=0.001;
%r1=r(1:end-1); r2=r(2:end);
%FQ(end,:)=N*pi*(r2.^2-r1.^2);
%fq(end,:)=[0 N*pi*(rm.^2-r1.^2)]+[N*pi*(r2.^2-rm.^2) 0];
%% Test steady state solution
phi0=Q/(2*pi*T)*besselk(0,r/lambda);           % Glee

phi1=fdm2c(            r,y,[0;T/B(2)],c,[Inf;Inf]    ,IBOUND,IH,FQ,'axial'); 
phi2=fdm2NodeCentered (r,y,[1e-6;T/B(2)],[B(1)/c;1e6],ibound,ih,fq,0);
phi3=fdm2NodeCentered (r,y,[1e-6;T/B(2)],[B(1)/c;1e6],ibound,ih,fq,1');
phi4=fdm2NodeCentered (r,y,[1e-6;T/B(2)],[B(1)/c;1e6],ibound,ih,fq,2);

% fh=NaN(  Ny,Nr+1); fh(1,:)=0;
% fq=zeros(Ny,Nr+1); fq(2,1)=Q;

% phi2=fdm2NodeCentered (r,y(2:end),2*T/B(2),b/c,ibound,ih,fq,0);
% phi3=fdm2NodeCentered (r,y(2:end),2*T/B(2),b/c,ibound,ih,fq,1');
% phi4=fdm2NodeCentered (r,y(2:end),2*T/B(2),b/c,ibound,ih,fq,2);

figure; hold on; grid on; xlabel('r [m]'); ylabel('head [m]'); 
set(gca,'xscale','log');

title('Dupuit, Blue:analytic, Red:numeric');

q=Q/(2*pi*T);

%% Six times the bessel function
plot(logspace(-6,1,71),besselk(0,logspace(-6,1,71)),'c'); % Bessel function directly
plot(r /lambda, phi0/q,'k');  % analytic
plot(rm/lambda, phi1/q,'bx'); % fdm cell centered
plot(r /lambda ,phi2/q,'r');  % mesh centered
plot(r /lambda ,phi3/q,'gs'); % mesh centered
plot(r /lambda ,phi4/q,'mp'); % mesh centered

% plot(r /lambda, phi0,'k');  % analytic
% plot(rm/lambda, phi1,'bx'); % fdm cell centered
% plot(r /lambda ,phi2,'r');  % mesh centered
% plot(r /lambda ,phi3,'gs'); % mesh centered
% plot(r /lambda ,phi4,'mp'); % mesh centered

% Seems to be perfectly correct, however

set(gca,'yscale','log');

% shows that values near zero tend to be inaccurate.
% This can be solved by taking many more points along the r-axis.
% the model seems sensitive to this for small drawdowns.
% take the number of points beyond 2000 and the entire line
% will be accurately computed by the numeric model.