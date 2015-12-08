%% mf_analyze script for Check balance

% Samani

% TO 190112

close all
clear variables

load underneath
load name
load(basename);

%% Position of observation points
Obs=[1.83 0 3.2;
     1.83 0 0.8];

Idx = xyzindex(Obs,gr);
nObs = numel(Idx);

%% % Reading simulation data

B=readBud([basename,'.bgt']);
B=mf_Psi(B); Prange=ContourRange(B,50,'Psi');

H=readDat([basename,'.hds']); Hrange=ContourRange(H,50);

phiObs=NaN(nObs,length(H));
for it=1:length(H)
    phiObs(:,it)= H(it).values(Idx);
end

%%
figure; hold on; set(gca,'xscale','log'); grid on;

%% Time drawdown for the two observation wells
title('time drawdown for the two observation wells');
xlabel('t [s]'); ylabel('ddn [m]');
plot([H.totim],phiObs);
legend('obs1','obs2');
