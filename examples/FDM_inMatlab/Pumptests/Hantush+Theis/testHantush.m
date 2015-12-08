%% Testing the standard Hantush function for transient flow in a confined aquifer
% TO 12-1-14

%%  Reproduces the tables of W(u,r/L) in Kruseman & De Ridder (1970)

% Annex IV block 1
rho=[0 0.005 [1 2 3 4 6 8]*10.^-2];
u=[1 2 4 6 8]' *10.^-(6:-1:1); u=[0 ; u(:)];

fprintf('\n');
fprintf('Kruseman & De Ridder (1970), Annex IV, block1, Hantush''s W(u,r/L)\n');
fprintf('First line [NaN r/L]\n');
fprintf('Following lines [u W(u,r/L]]\n');
Block1=Wh(u,rho);
display([[NaN; u] [rho; Block1]]);

% Annex IV block 2
rho=[0 [1 2 3 4 6 8]*10.^-1]; 
u=  [1 2 4 6 8]'*10.^-(4:-1:0); u=[0; u(u<=2)];

fprintf('\n');
fprintf('Kruseman & De Ridder (1970), Annex IV, block2, Hantush''s W(u,r/L)\n');
fprintf('First line [NaN r/L]\n');
fprintf('Following lines [u W(u,r/L]]\n');
Block2=Wh(u,rho);
display([[NaN; u] [rho; Block2]]);

% Annex IV block 3
rho=[1 2 3 4 6];
u  =[1 2 4 6 8]' *10.^-(2:-1:0); u=[0; u(u<=4)];

fprintf('\n');
fprintf('Kruseman & De Ridder (1970), Annex IV, block3, Hantush''s W(u,r/L)\n');
fprintf('First line [NaN r/L]\n');
fprintf('Following lines [u W(u,r/L]]\n');
Block3=Wh(u,rho);
display([[NaN; u] [rho; Block3]]);



%% Produces Hantush type curves through function Wh(u,rho)

u=logspace(-6,1,10*(1-(-6)))';
rho=[1 2 3 4 5 7]'*10.^(-1:0); rho=[0; rho(:)];

w=Wh(u,rho);

figure; hold on; grid on;
title('Hantush type curves'); xlabel('1/u'); ylabel('Wh(u,rho)');
set(gca,'xscale','log','yscale','log');

j=find(1./u>1e4,1,'last');
xLbl=1e5;
for i=2:length(rho)
    plot(1./u,w(:,i));
    h=text(xLbl,w(j,i),sprintf('r/L=%g',rho(i)),'BackgroundColor','W','fontsize',12);
    xLbl=xLbl/10; if xLbl<1e2, xLbl=1e5; end
end
%
i=floor(0.85*length(u));
plot(1./u,w(:,1),'k','linewidth',2);
text(0.5*(1./u(i)),w(i,1),'Theis','rotation',45,'fontsize',12,'BackgroundColor','w')
