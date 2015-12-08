%% mf_analyzeALL.m
% Example see USGS modpath Version 6 (2012), all examples named
% simulation1 through simulation5, see manual p36ff
% TO 130220

% TO 130221
clear variables; close all

load('name.mat') % get basename stored in file name.mat
load(basename);  % having retrieved baasename load the data in basename.mat
load underneath  % to get gr object

H=readDat([basename,'','.hds']); % read the unformatted head file

hrange=ContourRange(H,100); % get contour intervals
grey=[0.8 0.8 0.8]; % background color from figure

fsz=10; % fontsize

% for all layers make contour plot
ax = NaN(gr.Nlay,1);

%theseLayers = [1 5];
if ~exist('theseLayers','var'), theseLayers = [1 5]; end
    
for iLay=theseLayers
    figure; ax(iLay)=axes('nextplot','add','fontsize',fsz);
    set(ax(iLay),'clim',hrange(([1 end]),...
        'xlim',gr.xGr([1 end]),'xgrid','on',...
        'ylim',gr.yGr([end 1]),'ygrid','on',...
        'zlim',gr.zGr([end 1]),'zgrid','on'); 
    
    xlabel(ax(iLay),'x (ft)');
    ylabel(ax(iLay),'y (ft)');
    title( ax(iLay), sprintf('Modpath6 %s, head in layer %d',basename,iLay))
    
    hrange=ContourRange(H,500,iLay); % get contour intervals

    [c,hdl]=contour(ax(iLay),gr.xm,gr.ym,H(end).values(:,:,iLay),hrange,'b'); % contour layer iLay
    clabel(c,hdl); % plot contour labels

    % plot wells in this layer
    Iw=find([well.iLay]==iLay);
    for iw=1:length(Iw)
        well(iw).plotXY('or','parent',ax(iLay));
    end

    gr.plotBCN('RIV',RIV,3,1,'r','lineWidth',2);
    well.plotXY(ax(iLay));        
end

%% Use zone budget to get budget overview

B=readBud6([basename, '.bgt']); % read the Budget file
zonebudget(B);
