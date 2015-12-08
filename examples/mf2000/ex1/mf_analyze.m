%% Ex1 -- Analyzing output of the model
% Visualization of example EX1 of modflow 2000
%
% TO 091011 091129

%% Load model name and the data generated and saved by mf_adapt
load('name.mat') % get basename stored in file name.mat
load(basename);  % having retrieved basename, load the data in [basename '.mat']
load underneath  % load possible non-standard info from mf_adapt

fsz=10;          % set fontsize for plots

%% Read computed head data

H=readDat([basename,'.hds']);  % read the file

hrange = ContourRange(H,50);   % determine suitable contour elevations

%% Perpare plots

figure('name','example1','position',screenPos(1)); % new figure

P=[0 0 0.4 0.4]; % size of subaxis

% create axes at their desired position
ax(1)=axes('position', [0.05 0.55 0 0]+P,'nextplot','add','fontsize',fsz);
ax(2)=axes('position', [0.55 0.55 0 0]+P,'nextplot','add','fontsize',fsz);
ax(3)=axes('position', [0.05 0.05 0 0]+P,'nextplot','add','fontsize',fsz);
ax(4)=axes('position', [0.55 0.05 0 0]+P,'nextplot','add','fontsize',fsz,...
        'xColor',grey,'yColor',grey,'color','k','ticklength',[0 0]);

% for all layers make contour plot
for iLay=1:gr.Nlay
    set(gca,'clim',hrange([1 end]));
    
    xlabel(ax(iLay),'x [1000ft]');
    ylabel(ax(iLay),'y [1000ft]');
    title( ax(iLay), sprintf('head in layer %d',iLay))

    [c,hdl] = contourf(ax(iLay),gr.xm,gr.ym,H.values(:,:,iLay),hrange,'edgecolor','none'); % contour layer iLay
    clabel(c,hdl); % plot contour labels

    Iw=find([well.iLay]==iLay);
    well.plotXY(ax(iLay),'ro');

    if iLay == 1, plot(ax(iLay),drn(:,1),drn(:,2),'color','g','linewidth',2); end
    
end

%% Annotate plot with text
s={'Example from MODFLOW 2000 manual, p89';''; 'Computed in mfLab';''; 'TO 2012-04-07'};
text(0.5,0.5,s,'HorizontalAlignment','center','color','yellow','fontsize',fsz);

%% Use zonebudget to get budget overview

B=readBud([basename, '.bgt']); % read the Budget file

Zonearray=gr.const(1:gr.Nlay);  % first create a zone array (here onze zone per layer)
zonebudget(B,Zonearray,1); % budget of layer 1
zonebudget(B,Zonearray,2); % same for layer 2
zonebudget(B,Zonearray,3); % same for layer 3
zonebudget(B,Zonearray,[1 2 3]); % for all 3 layers totaled
