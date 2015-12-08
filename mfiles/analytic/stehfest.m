function [ft,sumKn,sumKnOvern]=stehfest(Fs,t,N)
%STEHFEST numerical back transformation from Laplace space according to Stehfest
%
% Example:
%    stehfest();  % is selfTest (uses Laplace transform of Theis' drawdown)
%    [ft,sumKn,sumKnOvern]=stehfest(Fs,t,N)
%     ft                 =stehfest(Fs,t);
%
%   Fs is the function in laplace space which accepts argument(s)
%
%   See also hantushn
%
%  TO 120120 (See Lee (1999))

if nargin==0
    selfTest();
    return;
end

if nargin<3, N=14; end % seems optimal for theis

ft=0;
sumKn=0;
sumKnOvern=0;
for n=1:N
    Kn=0;
    for k=floor((n+1)/2):min(n,N/2)
        Kn=Kn+k^(N/2)*factorial(2*k)./...
            (factorial(N/2-k)*factorial(k)*factorial(k-1)*factorial(n-k)*factorial(2*k-n));
    end
    Kn=(-1)^(n+N/2)*Kn;
    ft=ft+Kn*Fs(n*log(2)/t);  % n*log(2)/t is the Laplace parameter
    sumKn=sumKn+Kn;
    sumKnOvern=sumKnOvern+Kn/n;
end
ft=ft*log(2)/t;

end

function selfTest()
%SELFTEST tests stehfest when stehfest is calle without arguments
%
% The Laplace transform of the Theis solution is
% L(expint(kappa*t) = 2*K0(sqrt(kappa*s))/s
%
% TO 120120

    t=1;
    %kappa=r^2*S/4*kD;  -> sqrt(kappa)=sqrt(r^2S/kD)/2
    kappa=0.000001;
    u=kappa/t; % = r^2S/(4kDt)

    for N=2:2:40
        [W, sumKn, sumKnOvern]=stehfest(@(s) 2*besselk(0,2*sqrt(kappa*s))/s,t,N);

        fprintf('N=%2d: u=%12.7g, expint=%12.7g, W=%12.7g, rerr=%14.7f, aerr=%12g, sumKn=%12g, sumKnOvern=%12.7g\n',...
            N,u,expint(u),W,expint(u)/W-1,expint(u)-W,sumKn,sumKnOvern);
    end

end