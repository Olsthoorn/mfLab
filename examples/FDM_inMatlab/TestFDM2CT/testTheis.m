% This mfile runs a test of the finite element model fdm2t implemented
% directly in Matlab. The test is applied on axially symmetric flow.
% First the steady state variant, then the transient version
% Semi-confined flow (Glee and Hantush)
% TO 120110

close all
clear

%% Paremetrs
T=1000; c=1000; b=1; kh=T/b; kv=1e3; S =0.1; Q=2400;

lambda=sqrt(T*c);

r1=1e-6*lambda;
r2=1e2 *lambda;
n=50;

r=logspace(log10(r1),log10(r2),n);
y=[-b 0];
[r,y,rm,ym,dr,dy,Nr,Ny]=modelsize(r,y);

B=dy;

%% Numerical boundary conditions
IBOUND=ones(Ny  ,Nr  );    ibound=ones(Ny+1,Nr+1);

IH=zeros(Ny,Nr);           ih=zeros  (Ny+1,Nr+1);
FQ=zeros(Ny,Nr);           fq=zeros(Ny+1,Nr+1);
FQ(end,1)=Q;               fq(end-1:end,1)=Q/2;

%% Test steady state solution

t=[0 logspace(-2,2,41)];

phi0=Q/(4*pi*T)*expint(S./(4*T*t(:))*r.^2);           % Glee
phi1=fdm2ct(        r,y,t,kh,[],kv,S,IBOUND,IH,FQ,'axial'); 
phi2=fdm2tNC (r,y,t,kh,kv,S,ibound,ih,fq,0);
phi3=fdm2tNC (r,y,t,kh,kv,S,ibound,ih,fq,1');
phi4=fdm2tNC (r,y,t,kh,kv,S,ibound,ih,fq,2);

% ih=zeros(Ny,Nr+1);
% fq=zeros(Ny,Nr+1); fq(2,1)=Q;

% phi2=fdm2NC (r,y(2:end),2*T/B(2),b/c,ibound,ih,fq,0);
% phi3=fdm2NC (r,y(2:end),2*T/B(2),b/c,ibound,ih,fq,1');
% phi4=fdm2NC (r,y(2:end),2*T/B(2),b/c,ibound,ih,fq,2);

phi1=permute(phi1,[3,2,1]);
phi2=permute(phi2,[3,2,1]);
phi3=permute(phi3,[3,2,1]);
phi4=permute(phi4,[3,2,1]);

figure; hold on; grid on; xlabel('r [m]'); ylabel('head [m]'); 
set(gca,'xscale','log');

title('Dupuit, Blue:analytic, Red:numeric');

q=Q/(4*pi*T);

%% Six times the bessel function
%plot(logspace(-6,1,71),besselk(0,logspace(-6,1,71)),'c'); % Bessel function directly
plot(r , phi0/q,'k');  % analytic
plot(rm, phi1(:,:,end)/q,'bx'); % fdm cell centered
plot(r  ,phi2(:,:,end)/q,'r');  % mesh centered
plot(r  ,phi3(:,:,end)/q,'gs'); % mesh centered
plot(r  ,phi4(:,:,end)/q,'mp'); % mesh centered

% plot(r , phi0,'k','linewidth',2);  % analytic
% plot(rm, phi1(:,:,end),'bx'); % fdm cell centered
% plot(r  ,phi2(:,:,end),'r');  % mesh centered
% plot(r  ,phi3(:,:,end),'gs'); % mesh centered
% plot(r  ,phi4(:,:,end),'mp'); % mesh centered

% Seems to be perfectly correct, however

set(gca,'yscale','log');

% shows that values near zero tend to be inaccurate.
% This can be solved by taking many more points along the r-axis.
% the model seems sensitive to this for small drawdowns.
% take the number of points beyond 2000 and the entire line
% will be accurately computed by the numeric model.

%% Type curves

u =(S./(4*T*t(:)))*(r.^2);
WM1=phi4(:,:,end)/q;
WM2=phi0/q;

figure; hold on; 
for i=1:size(u,2)
    plot(1./u(:,i),WM1(:,i),'b');
    plot(1./u(:,i),WM2(:,i),'r');
    set(gca,'xscale','log','yscale','log');
end
xlabel('1/u'); ylabel('W(u)');