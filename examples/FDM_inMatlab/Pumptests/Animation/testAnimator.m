%TESTANIMATOR tests animator of on-screen pumping test analysis
%
% How to do a pumping test by shifting data on the screen
%
% See also: animator_loglog
%
% TO 121201

clear variables;
close all;

figure; hold on;

T=1000; c=100; Sa =0.00001; Sy=0.1; Q=2400; Alpha=0.5; c=Sy/Alpha;

%% Data from pumping test Vennebulten (Kruseman & De Ridder, 1970)

[vnams,vals]=getExcelData('MyData.xls','Vennebulten','Hor');


Qd=vals(:,strmatchi('Q',vnams));
td=vals(:,strmatchi('t',vnams));
sd=vals(:,strmatchi('ddn90/16',vnams));
rd=90;

Qd=Qd(~isnan(sd));
td=td(~isnan(sd))/(24*60);
sd=sd(~isnan(sd));

Q=Qd(1);


%% Show it

% Notice the ButtonDownFcn, the animator for log log scale and the
% erasemode

%% any object you want to be selectable to drag the group of objects gets the same properties

%% My data
plot(td,sd,'ro','tag','myData','ButtonDownFcn','animator_loglog start','EraseMode','xor');

%% A Match point initially on 1,1 (useful for pumping tests)
xMP=1;
yMP=1;
plot(xMP,yMP,'ko','tag','myData','ButtonDownFcn','animator_loglog start','EraseMode','xor');

set(gca,'xscale','log','yscale','log');
xlabel('1/u'); ylabel('s');
grid on


%% Now you can pick up the data on the figure with the mose and move it