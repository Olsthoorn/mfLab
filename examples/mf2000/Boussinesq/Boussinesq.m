%Boussenesq analytical Steward WRR 2007
% Analytical solution Boussinesq, Polibarinova Kochina (Steward, 2007)

slope=-1/2;  % tan(theta)
Q=5;         % m2/d
k=0.2;       % m/d

H0=-Q/(k*slope); %normal depth

% eta = H/H0, eta1=H1/H0
eta1A=2.0; etaA = logspace(0,log10(4),50);
eta1B=0.5; etaB = linspace(0.001,0.999,100);
sA=H0/(-slope)*(log((etaA-1)/(eta1A-1))+etaA-eta1A);
sB=H0/(-slope)*(log((etaB-1)/(eta1B-1))+etaB-eta1B);

xlabel('x-x0 [m]'); ylabel('phi'); title('Boussinesq analytical');
if 0
plot(sA,H0*etaA+sA*slope,'b',sB,H0*etaB+sB*slope,'r',...
    sA,sA*slope,'k',sB,sB*slope,'k',...
    sA,sA*slope+H0,'g',sB,sB*slope+H0,'g');
else % dimensionless
plot(sA/H0*(-slope),etaA,'b',sB/H0*(-slope),etaB,'r',...
     sA/H0*(-slope),ones(size(sA)),'g',...
     sB/H0*(-slope),ones(size(sB)),'g');
xlabel('(x-x_1) tan(-\theta)/H_0');
ylabel('\eta=H/H_0');
title('Analytical Boussinesq solution, for H=2H_0 and H=0.5H_0 at x=x_1')
end