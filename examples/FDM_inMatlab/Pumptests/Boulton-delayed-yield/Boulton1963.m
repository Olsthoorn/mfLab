% Boulton delaye yield
% TO 120112

if ~exist('VENNEBULTEN','var')
    close all
    clear

%% Paremetrs
    T=1000; c=100; Sa =0.00001; Sy=0.1; Q=2400;
end

%% Grid
b=1; 
lambda=sqrt(T*c);
r1=1e-6*lambda;
r2=1e3 *lambda;
r=logspace(log10(r1),log10(r2),10*(log10(r2)-log10(r1)));
y=[-b 0 1];
[r,y,rm,ym,dr,dy,Nr,Ny]=modelsize(r,y);

B=dy;

%% The model
kh=[0  ;T/B(2)];
kv=[B(1)/c; Inf];
S =[Sy/B(1); Sa/B(2)];
alpha =1/(c*Sa);

%% Numerical boundary conditions
IBOUND=ones(Ny,Nr);
IH=zeros(Ny,Nr);     IH(end,end)=0;
FQ=zeros(Ny,Nr);
FQ(end,1)=Q;

%% Test steady state solution
t1=1e-5;
t2=1e3;
t=logspace(log10(t1),log10(t2),10*(log10(t2)-log10(t1)))';

I=find(t>=1e-4); % remove first times

phia=Q/(4*pi*T)*expint(Sa./(4*T*t(:))*r.^2);
phiy=Q/(4*pi*T)*expint(Sy./(4*T*t(:))*r.^2);

IH(end,:)=Q/(4*pi*T)*expint(Sa./(4*T*t(1))*rm.^2);

phiB=fdm2t(r,y,t,kh,kv,S,IBOUND,IH,FQ,'axial'); phiB=permute(phiB,[3,2,1]);

figure; hold on; grid on; fontsize=15;
xlabel('1/u_a = t/r^2 * (4T/S_y)','fontsize',fontsize);
ylabel('s/(4\pi T)','fontsize',fontsize); 
title('Boulton type curves','fontsize',fontsize);

ylim=[1e-2 1e2]; xlim=[1e-5 1e4];
set(gca,'xscale','log','yscale','log','ylim',ylim,'xlim',xlim,'fontsize',fontsize);

q=Q/(4*pi*T);

gamma=1+Sy/Sa; cutoff = 3 /gamma;

ua=Sy./(4*T*t)*r .^2; ua(ua>3*gamma)=NaN;
uy=Sy./(4*T*t)*r .^2; uy(uy>3*gamma)=NaN;
uB=Sy./(4*T*t)*rm.^2; uB(uB>3*gamma)=NaN;
plot(1./ua(I,:), phia(I,:)/q,'k');  % analytic
plot(1./uy(I,:), phiy(I,:)/q,'k');
plot(1./uB(I,:), phiB(I,:,end)/q,'b');

BB=sqrt(rm/lambda); BBstr=sprintf('%10.1f\n',BB);

i=round(length(I)/2);
for j=1:size(uB,2)
    h=text(Sy./uB(I(i),j)',phiB(I(i),j,end)'/q,sprintf('%.2f',BB(j)),'clipping','on');
end
