%% This file use used to visualize the results produced by MODFLOW after
% mf_setup finised and the model has finished normally
% Both the analytic and the numerical model show the idial drawdown
% is about 0.7 m, while the drawdown in reality is over 3.5 m. Therefore
% quite some development of this HDDW needs still to be done.
%
%  TO 100610 130222

close all
load HDDW        % model data saved in mf_adapt
load underneath  % stuff we wanted to remember from mf_adapt

% Read the heads produced by MODFLOW
H=readDat([basename '.HDS']);

%% Show results
figure; hold on;
hrange = ContourRange(H,50);
contourf(gr.xc,gr.yc,H(end).values(:,:,well(end).iLay(1)),'k');  % hrange,'k');
set(gca,'xlim',[0 50],'ylim',[-25 25]);
title('xy plane'); xlabel('x [m]'); ylabel('y [m]');

%%
figure;
contourf(gr.xc,gr.zc ,XS(H(1).values(hit(gr.yGr,0),:,:)));
title('zx plane'); xlabel('x'); ylabel('z');
set(gca,'xlim',[0 20],'ylim',[-20 0]);

%%
figure;
contourf(gr.yc,gr.zc,YS(H(1).values(:,1,:)));
title('yz plane'); xlabel('y'); ylabel('z');
set(gca,'xlim',[-20 20],'ylim',[-20 0]);

%%
figure;
idx = well(1).idx;

plot(gr.XM(idx),H(end).values(idx),'b','linewidth',3);
title('Computed head in ideal 30 m long HDDW for 19 m3/h extraction');
xlabel('coordinate along HDDW from center');
ylabel('drawdown [m]');

%% Check HDDW discharge
B=readBud([basename '.BGT']);

Lbl = 'WELLS';

for it=1:numel(H)
    fprintf('Q_HDDW(it=%d,t=%g d) = %g m3/d. ',it,H(it).totim,sum(B(it).term{strmatchi(Lbl,B(it).label)}(idx)));
    fprintf('Havg in HDDW = %g m\n',mean(H(it).values(idx)));
end

