%DTOPBALANCE Processing water balance for spreadsheet
%
% Works for examples in mflab/examples/fm2005/Dutchtop
% it may be specific to the Dutch situation.
%
% Might also be useful as a startup for other examples
%
% TO 100925

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

clear;  % close all don't close to compare with other figure

load name  % get basename
load Balance.mat

fprintf('Loading and processing water Balance\n');
fprintf('File was saved at %s\n',datestr(timestamp));

%% saving the waterbalance for all sections
fprintf('Labels in budget file used for this balance\n');
for i=1:length(labels), fprintf('%-20s',labels{i}); end; fprintf('\n');

fprintf('\nMatrices corresponding to water balance terms\n');
for i=1:length(legtxt)
    fprintf('%-10s = %-10s\n',legtxt{i},labtxt{i});
end

fprintf('\nThe size of each water balance matrix is %d parcels by  %d days.\n',size(QWEL));

IN ={max(QWEL,0),max(QSTO,0),max(QGHB,0),max(QDRN,0),max(QRCH,0)};
OUT={min(QWEL,0),min(QSTO,0),min(QGHB,0),min(QDRN,0),min(QRCH,0)};

IN3 =zeros([size(QWEL) length(IN)]);
OUT3=zeros([size(QWEL) length(IN)]);

for i=1:length( IN), IN3( :,:,i)=IN{ i}; end
for i=1:length(OUT), OUT3(:,:,i)=OUT{i}; end

%% plot overall water balance

figure; hold on;
grid on; xlabel(sprintf('%s - %s',datestr(tne(1,1),'yyyy'),datestr(tne(end,1),'yyyy'))); ylabel('m/d');

title(sprintf('running water balance for %s weighted over all parcels, total area=%.0f ha\n',...
    basename,sum([P.area])/1e4));

%% Change to weighted value over all sections so that total becomes in m/d

areaweight=repmat(([P.area]'*ones(1,size(QWEL,2)))/sum([P.area]),[1,1,length(IN)]);

in3 = areaweight.*IN3 ;  % in m/d
out3= areaweight.*OUT3;  % in m/d

%% plot running water balance in m/d
area(tne(:,1),sum(permute(in3 ,[3,2,1]),3)');  % summing over all sections (first dimension)
area(tne(:,1),sum(permute(out3,[3,2,1]),3)'); % summing over all sections (first dimension)
legend(legtxt);
datetick('x',4);

%% Save the water balance values in to csv files, and timestamp the files

T=datevec(tne(:,1));

dlmwrite([basename '-total-in--UPSPG-STO-DITCH-ROFF-EVT ' datestr(timestamp,'dd-mmm-yyyy HH:MM')],[T(:,1:3) squeeze(sum(in3, 1))],'delimiter','\t');
dlmwrite([basename '-total-out-DNSPG-STO-DITCH-RIN--RCH ' datestr(timestamp,'dd-mmm-yyyy HH:MM')],[T(:,1:3) squeeze(sum(out3,1))],'delimiter','\t');

fprintf('Done\n');

