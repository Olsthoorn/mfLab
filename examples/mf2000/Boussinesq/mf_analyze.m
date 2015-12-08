%% Analyzing output of the model
% TO 100415 100916
 
load('name.mat') % loads name.mat, which only contains the variable "basename"
load(basename);  % having retrieved the baasename value, we can load

%% load the unformatted head file
H=readDat([basename,'','.hds']);  % use readDAT to read the heads  output file
H=maskHC(H,1000);                 % throws out HNOFLO and HDRY and inactive cells


%% Plot bottom of aquifer and start heads
figure; hold on
xlabel('x [m]'); ylabel('elevation [m]'); grid on;
title('Aquifer bottom and head elevation on hill slope');

plot(gr.xm,gr.Z(:,:,end),'color','k','linewidth',3); % plot bottom of aquifer thick
plot(gr.xm,STRTHD(:,:,end),'r');                  % plot start heads in read

%% Plot heads at various times as head snapshots in gread
for it=1:length(H)
    for iL=1:size(H(1).values,3)
        plot(gr.xm,H(it).values(:,:,iL),'g');
    end
end

plot(gr.xm,H(end).values(:,:,iL),'r');

%% Read unformatted budget file and mask noflow cells if they exist

B=readBud([basename,'','.bgt']);
zonebudget(B)
