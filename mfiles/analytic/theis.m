function w=theis(u)
%THEIS Drawdown according to Theis computed as hantush(u,0) + selftest compares integration with expint.
%
% Example:
%    w=theis(u);
%
% Theis well function by integration
% TO 010409 101211

if nargin<1, selftest; return; end

w=hantush(u,0);

function w=selftest
%SELFTEST tests itself
%
u = logspace(-6,1,71);
w = hantush(u,0);

figure; loglog(1./u,expint(u),'r','linewidth',1); hold on

plot(  1./u,        w,'o'); xlabel('1/u'); ylabel('W');

grid on; set(gca,'xlim',[1e-1 1e6],'ylim',[1e-5 100],'color','none');
legend('expint','hantush(u,0)');
title('Theis'' well function');
