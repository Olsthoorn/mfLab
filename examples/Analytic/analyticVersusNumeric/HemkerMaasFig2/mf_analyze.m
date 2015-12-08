%% Analyzing output of the model
% 131005
 
load('name.mat') % get basename stored in file name.mat
load(basename);  % having retrieved baasename load the data in basename.mat
load underneath  % to get gr object

ir= hit(gr.xGr,r);  % get column for piezometer used by HM (1987_

%% Analytical n-layer drawdown
IL = kD(:,1)>1;
drawdown = hantushn(Q,r,t,St,Sf(2:end-1),c,kD(IL));
drawdown = squeeze(abs(drawdown))';

%% Heads (drawdown) from MODFLOE)
H=readDat([basename,'','.hds']); % read the unformatted head file

% Select only the aquifers
IL = find(HK(:,1,:)>1);

% select the computed drawdowns for piezometer at r and all times 
drawdownMdl = NaN(size(drawdown));  % MODLFOW
drawdownFDM = NaN(size(drawdown));  % fdm2ct
for it=1:numel(H)
    drawdownMdl(it,:) = -H(it).values(1,ir,IL);   % layers vary due to aquitards subdivision
    drawdownFDM(it,:) = -Phi([3 5],ir,it);        % always 5 layers
end

% print outcomes of analytic, MODFLOW and fdm2ct models
fprintf('%10.3g  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f\n',[t drawdown drawdownMdl drawdownFDM]');

%% Plots
figure; hold on; grid on;

plot(t,drawdown,'r','linewidth',2);      leg1 = {'Hemker/Maas (1987) L1', 'idem L2'};
plot(t,drawdownMdl,'bo-');               leg2 = {'MODFLOW L1','MODFLOW L2'};

%plot(t,drawdownFDM,'gx-'); leg3 = {'fdm2ct L1','fdm2ct L2'}; % from fdm2ct

plot(t,sTheis1 ,'m');                    leg4 = {'Theis L2'};
%plot(t,sTheis12(:,1),'c');              leg4 = {'Theis L2'}; % 'Theis L1+2'};

plot(t,sHantush1,'mo-');                 leg5 = {'Hantush L2'}; %,'Hantush L1+2'};
%plot(t,sHantush12,'co-');               leg5 = {'Hantush L2','Hantush L1+2'};

set(gca,'xscale','log');

legend(leg1{:},leg2{:},leg4{:},leg5{:},2); % leg3{:}, leg5{:}

set(gca,'ylim',[0 6]);

% Plot vertial line at characteristic time of this problem (Sc/4=1.6e-3*100/4)
plot([0.04 0.04],get(gca,'ylim'),'r--');

xlabel('time [d]'); ylabel('drawdown [m]');
title(sprintf('Hemker/Maas (1987) vs MODFLOW, %d layers per aquitard',layersPerAquitard));


%% Plot on double log paper directly on the original image by Hemker and Maas, 1978 (fig 2)

% Hand-digitized pixel coordinats of figure in image HemkerMaas87-Fig2.png.
pLL=[   -0.289565     243.552 ];  % pixels, whole picture
pUR=[     241.71      0.289565];  % pixsls, whole picture

fLL=[      20.333      111.399];  % pixels axis of graph
fUR=[     226.138      7.44435];  % pixels axis of graph

LL = [0.1 0.01]; LfLL = log10(LL);  % graph coordinates in graph units and their log10
UR = [1e5 10];   LfUR = log10(UR);  % same
  
% Get the log10 coordinates of the extends of the whole image to allow
% putting it on correct axis
Lx = interp1([fLL(1) fUR(1)],[LfLL(1) LfUR(1)],[pLL(1),pUR(1)],'linear','extrap');
Ly = interp1([fLL(2) fUR(2)],[LfLL(2) LfUR(2)],[pLL(2),pUR(2)],'linear','extrap');

%% now plot directly on the image
figure;
A = imread('HemkerMaas87-Fig2.png');                  % get image
image(Lx,Ly(end:-1:1),A); set(gca,'ydir','normal')    % put it on axis on coordinates

hold on;

% find computed drawdowns that fit within the graph on the image
% analytic drawdowns, and plot them
I1 = drawdown(:,1)>0.5e-2; plot(log10(t(I1))+4,log10(drawdown(I1,1)),'r','lineWidth',2);
I2 = drawdown(:,2)>0.5e-2; plot(log10(t(I2))+4,log10(drawdown(I2,2)),'r','lineWidth',2);
% modflow drawdowns, and plot them
I1 = drawdownMdl(:,1)>0.5e-2; plot(log10(t(I1))+4,log10(drawdownMdl(I1,1)),'bo-');
I2 = drawdownMdl(:,2)>0.5e-2; plot(log10(t(I2))+4,log10(drawdownMdl(I2,2)),'bo-');

