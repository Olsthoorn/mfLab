% grafiek Gamma
% Gamma = (b/lambda) coth (b/lambda) - 1

b_L = logspace(-2,2,41);
y = b_L.*coth(b_L)-1;
ylim = (1/3)*b_L.^2;
figure; hold on;
plot(b_L,y)
title('y = b/\lambda coth( b/\lambda ) - 1')
xlabel('b / \lambda');
ylabel('y');
set(gca,'xScale','log','xGrid','on','yGrid','on');