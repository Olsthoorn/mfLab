function TypeCurvesHantush
%TYPECURVESHANTUSH produces Hantush type curves through function Wh(u,rho)
%
%
% TO 120114

%% Suitable ranges
u=logspace(-6,1,10*(1-(-6)))';
rho=[1 2 3 4 5 7]'*10.^(-2:0); rho=[0; rho(:)];

w=Wh(u,rho); % compute Hantush values for all combinations

figure; hold on; grid on;
title('Theis and Hantush type curves','fontsize',16);
xlabel('1/u','fontsize',16);
ylabel('Wh(u,rho)','fontsize',16);
set(gca,'fontsize',14,'xscale','log','yscale','log',...
    'xGrid','on','ygrid','on',...
    'xminortick','on','yminortick','on');

j=find(1./u>1e4,1,'last');
xLbl=1e5;
for i=2:length(rho)
    plot(1./u,w(:,i));
    h=text(xLbl,w(j,i),['r/\lambda =' sprintf('%.2f',rho(i))],'BackgroundColor','W','fontsize',12);
    xLbl=xLbl/10; if xLbl<1e2, xLbl=1e5; end
end

i=floor(0.85*length(u));
plot(1./u,w(:,1),'k','linewidth',2);
text(0.5*(1./u(i)),w(i,1),'Theis','rotation',45,'fontsize',12,'BackgroundColor','w')

