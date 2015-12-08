function [Wh1,u,rho]=Wh1(u,rho)
%WH1 computes Hantush's well function (Maas C, Veling, E, (2010), Stromingen, 16(2010)59-69
%
% Example:
%    [W,u,rho]=Wh1(u,rho)
%
% INPUT:
%     u   = r^2S/(4kD t) as a vector. It will be ordered vertically
%     rho = r/lambda as a vector. It will be ordered horizontally
% OUTPUT
%     W   = well function, time vertially distance horizontally
%     u   = same size and orientation as W
%     rho = same size and orientation as W
%
% In their paper Maas and Veling use rho and tau, where rho is r/lambda and
% tau = ln (2/rho * t/(c*S). This tau can be converted back to traditional
% parameteres rho and u:
% tau = ln (rho/2 * 4 kD t/(r^2 S))
%     = ln (rho/2 / u)
%
% In this function, we assume that r is oriented horizontally and
% t vertically. To transpose, use use Wh(u,rho)'
%
% See also: Wh, hantush, hantushn
%
% TO 13042

%% Allow vector input to cover any combination of a list of u and of rho (r/L) values

if nargin==0, selfTest(); return; end

if nargin<2, rho=0; end

rho = rho(:)';
u   = u(:);

tau = log(bsxfun(@rdivide,rho(:)/2,u'))';
tau(tau>100)=100;

[Nt,Nr] = size(tau);

Rho   = ones(Nt,1) * rho;
h_inf = besselk(0,Rho);
w     = (expint(Rho)-h_inf)./(expint(Rho)-expint(Rho/2));
I     = h_inf-w.*expint(Rho/2.*exp(abs(tau))) + (w-1).*expint(Rho.*cosh(tau));
Wh1   = h_inf +sign(tau).*I;

end

function selfTest

rho = [2 3 5 7.5 10]' * [0.01 0.1 1];
rho = rho(:);

u   = (2:10)' * logspace(-5,0,6);
u   = u(:);

[w,u,rho] = Wh1(u,rho);

figure; xlabel('1/u'); ylabel('Wh(u,r/lambda)');
title('Hanush''s Wh according to Maas and Veling (Strominge 16(2010)59-69');
set(gca,'nextPlot','add','xScale','log','yScale','log','xGrid','on','yGrid','on')

for ir = numel(rho):-1:1
    leg{ir} = sprintf('r/L = %g',rho(ir));
    plot(1./u,w(:,ir),mf_color(ir));
end
plot(1./u,expint(u),'ko');
leg=['Theis' leg];
legend(leg{end:-1:1},4);

end
