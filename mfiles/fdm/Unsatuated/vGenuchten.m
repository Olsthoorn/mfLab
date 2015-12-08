% Experiment Van Genuchten curves
%
% experiment with the shape of Van Genuchten curves to find a good
% parameter set for dedicated soils that can fulfill special functions like
% representing ponding.

% Voor een pond van 10 cm diep (kan 100 mm regen bergen) is n=3 geschikt)
% om deze hoeveelheid in 1 dag te kunnen infiltreren met gradient 1, dan
% moet Ks>=10 cm/d zijn.
close all; clear

dBodem = 10/2;
Ks     = 100;

h      = dBodem:.2:20;


% Van genuchten
n      = logspace(0,log10(20),10);
m      = 1-1./n;
L      = 0;
a      = 1/dBodem;
Se = @(n) (1+(a.*h).^n).^(-(1-1/n));
K  = @(m) Ks * Se(1./(1-m)).^L .* (1 - (1-Se(1./(1-m)).^(1./m)).^m).^2;

% Corey
psiAE = dBodem;
b = 0:2:10;
b = [0.1 0.3 1 3 10];
c = 2*b+3;
Sc = @(b) (h/psiAE).^(-1/b);
Kc = @(b) Ks * Sc(b).^(2*b+3);

figure;
defaults = {'nextplot','add','xGrid','on','yGrid','on','yScale','lin','ylim',[0 10],'yDir','reverse'};
ax1 = subplot(1,2,1,defaults{:});
ax2 = subplot(1,2,2,defaults{:},'xScale','log','xlim',[0.01 Ks]); 
title(ax1,'Saturation Se(h)');
title(ax2,'Conductivity K(h)');

xlabel(ax1,'Se [-]');
ylabel(ax1,'|h| [cm]');
ylabel(ax2,'|h| [cm]');
xlabel(ax2,'K [cm/d]');

leg{numel(n)}='fake';
for in = 1:numel(n)
    plot(ax1,Se(n(in)),h);
    m = 1-1./n(in);
    plot(ax2,K(m),h);
    leg{in} = sprintf('n=%g',n(in));
end
legend(ax1,leg{:});
legend(ax2,leg{:});

%%
figure;
ax3 = subplot(1,2,1,defaults{:});
ax4 = subplot(1,2,2,defaults{:},'xScale','log','xlim',[0.01 Ks]); 
title(ax3,'Saturation Corey & Campbell Sc(h)');
title(ax4,'Conductivity Corey & Campbell Kc(h)');

xlabel(ax3,'Sc [-]');
ylabel(ax3,'|h| [cm]');
ylabel(ax4,'|h| [cm]');
xlabel(ax4,'Kc [cm/d]');

legc{numel(b)}='fake';
for in = 1:numel(b)
    plot(ax3,Sc(b(in)),h);
    plot(ax4,Kc(b(in)),h);
    legc{in} = sprintf('b=%g',b(in));
end
legend(ax3,legc{:});
legend(ax4,legc{:});