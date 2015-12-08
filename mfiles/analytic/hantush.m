function w=hantush(u,rho)
%HANTUSH computes Hantush's well function computed by integration
%
% Example:
%    w=Hantush(u,r/Lambda)
%
% See also: hantushn hantushe Wh Wh1
%
% TO 010409

if nargin<2, w=selftest; return; end

if any(u)<=0;
    error('u must be >0');
end

u   = u(:);     % time (u) vertical
rho = rho(:)';  % distance (rho) horizotnal

onesr = ones(size(rho));

r2      = (rho/2).^2;
dksi    = 0.1;
ksi_end = 20;
ksi     = log(u);
w       = exp(-exp(ksi)*onesr-exp(-ksi)*r2).*dksi/2;  %initiation

while any(ksi<ksi_end)                        % Simpsons trapezium rule
    ksi = ksi+dksi;
    w   = w+exp(-exp(ksi)*onesr-exp(-ksi)*r2)*dksi;
end

function w=selftest
u=logspace(-6,1,71);
rho=[0.01 0.05 0.1 0.5 1 2 3 4 5 6];  % see Kruseman & De Ridder (1971)

w=hantush(u,rho);

figure;
n=0; leg={};
loglog(1./u,expint(u),'linewidth',3); n=n+1; leg{n}='Theis'; hold on
plot(  1./u,w); xlabel('1/u'); ylabel('W'); hold on

for i=1:length(rho)
    n=n+1; leg{n}=sprintf('r/L=%g',rho(i));
end
grid on; set(gca,'xlim',[1e-1 1e6],'ylim',[1e-5 100],'color','none');
legend(leg,4);

title('Hantush type curves self test');