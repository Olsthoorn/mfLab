function [H,HWanda,fCorr] = CJ13_8(R,lambda,rhoc,DT,t,d)
%CJ13_8_analyticHeatLoss --- computes analytic heatloss
%
% USAGE: 
%     H = CJ13_8_analyticHeatLoss(R,lambda,rhoc[,DT[,t[,d]]])
%     H= CJ13_8(R,lambda,rhoc)
% Out:
%   H [W/m] is heat loss from the surface of a cylinder of radius R of constant
%   temperature DT=1, in inifnite solid of heat conductivity lambda.
% In:
%   R = radius of pipe
%   lambda = heat conductance of porous medium (including water)
%   rhoc   = heat capacity    of porous medium (including water)
%   DT     = temperature change (default = 1)
%   t      = time, if omitted, t = R^2/kappa*tau,
%            where tau= logspace(-2,3,80) and kappa = lambda/rhoc
%   d      = depth of center of pipe below ground surface (d>R)
%            default = Inf, i.e. medium is infinite
%            If given, H is corrected for depth of pipe below ground
%            surface using H = H/(1+fCorr(t,2d))
%
% TO 121106
fsz = 16;

kappa = lambda/rhoc; % [J/s/m/K /[kg/m3]/[J/kg/K] = [m2/s]

% therefore, u has dimsnion [1/m] !!

j=100;  % span 10^-j to 10^j !! wide span is necessary for sufficient accuracy

% Trampizium rule for integration later on
u = logspace(-j,j,2000); du=diff(u); u = 0.5*(u(1:end-1)+u(2:end));

% Times to match figure 43 in Carslaw and Jaeger

if nargin<5
    tau = logspace(-2,3,80); t = R^2/kappa*tau;
else
    tau = kappa * t / R^2;
end

gamma = 0.57722;  % Euler's constant

if nargin>5 % than assume half-infinite porous medium instean of infinite porous medium
    r=2*d;
    if r<=R, error('d=%g must be > R=%g'); end
end

rho = r/R;

for it=length(tau):-1:1
    % complete solution
    %fCorrA(it) = CJp335_6A(rho,tau(it));  % not correct
    fCorr(it) = CJp335_6(rho,tau(it));
    f(it)     = CJp336_8(lambda,tau(it));
    
    
    % this gives f at the circumference
    %f( it) = 4/pi^2 * sum(arg1);
    
    if nargin>5
        fCorr(it) = CJp335_6(rho,tau(it));
    else
        fCorr(it)   = 0;
    end
    
    
    % solution for small times
    f1(it)= (1./sqrt(pi*tau(it))+(1/2)-sqrt(tau(it)/pi)/4+tau(it)/8);
    
    % solution for large times
    f2(it)= 2*(1./(log(2*tau(it))-2*gamma)-gamma./(log(4*tau(it))-2*gamma).^2);
    
    fWanda(it) = 1/log(1+sqrt(pi*tau(it)));
end

% remove instabilities/inaccuracies for extremely low tau
fCorr(  1:find(fCorr  <0,1,'last'))=0;  % instabilities

if nargin<4, DT=10; end

HWanda = 2*pi * lambda * DT * fWanda; 
H      =                 DT * f;
t      = R^2/kappa*tau;

figure; hold on; xlabel('time [s]'); ylabel('H [W/m]');
title(sprintf('Heat loss from pipe R=%g m and DT=%g C in infinite soil',R,DT));
plot(t,H); set(gca,'xscale','log','yscale','log','xgrid','on','ygrid','on');
set(gca,'fontSize',fsz)

% set up figure
figure;hold on;
xlabel('log_{10} \tau','fontsize',fsz);
ylabel('log_{10}(h R/ (\lambda T))','fontsize',fsz);
title('Fig 43 of Carslaw and Jaeger, p338','fontsize',fsz)
set(gca,'yscale','lin','xscale','lin','xgrid','on','ygrid','on',...
    'xlim',[-2 3],'ylim',[-1 1],'fontsize',fsz);


%% Analytical integral
% Plot the 3 lines
h=[];
h(1) = plot(log10(tau),log10(f/(2*pi*lambda)),'k');
h(2) = plot(log10(tau),log10(f1),'b');
h(3) = plot(log10(tau),log10(f2),'r');
h(4) = plot(log10(tau),log10(fWanda./(1-fCorr)),'m-');
h(5) = plot(log10(tau),log10(f),'k');

text(-1,-0.6,'\tau = \kappa t / R^2','fontsize',fsz,'backgroundColor','y');

legend(h,'CJ-13.8  (all   times)','CJ-13.9  (small  times)','CJ-13.10 (large times)','half-space','Inf-space');
fprintf('');


