%TESTSTEHFEST test of stehfest's numerical back transformation from Laplace space
%
% testStehfest is embedded as selfTest in stehfest
%
% see also: hantushn stehfest
%
% TO 120120

t=1;
%kappa=r^2*S/4*kD;  -> sqrt(kappa)=sqrt(r^2S/kD)/2
kappa=0.000001;
u=kappa/t; % = r^2S/(4kDt)

for N=2:2:40
%    The Laplace transform of the Theis solution is
%   L(expint(kappa*t) = 2*K0(sqrt(kappa*s))/s
    [W, sumKn, sumKnOvern]=stehfest(@(s) 2*besselk(0,2*sqrt(kappa*s))/s,t,N);

    fprintf('N=%2d: u=%12.7g, expint=%12.7g, W=%12.7g, rerr=%14.7f, aerr=%12g, sumKn=%12g, sumKnOvern=%12.7g\n',...
        N,u,expint(u),W,expint(u)/W-1,expint(u)-W,sumKn,sumKnOvern);
end