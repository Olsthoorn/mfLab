%% Analyzing output of the model
% TO 091011 091129 120413
 
load('name.mat') % get basename stored in file name.mat
load(basename);  % having retrieved baasename load the data in basename.mat
load underneath  % to get gr object

H=readDat([basename,'','.hds']); % read the unformatted head file

hrange = ContourRange(H,50); % get contour intervals

figure; % new figure
P=[0 0 0.4 0.4]; % size of subaxis

grey=get(gcf,'color'); % get background color from figure
fsz=10;

% create axes at their desired position
ax(1)=axes('position', [0.05 0.55 0 0]+P,'nextplot','add','fontsize',fsz);
ax(2)=axes('position', [0.55 0.55 0 0]+P,'nextplot','add','fontsize',fsz);
ax(3)=axes('position', [0.05 0.05 0 0]+P,'nextplot','add','fontsize',fsz);
ax(4)=axes('position', [0.55 0.05 0 0]+P,'nextplot','add','fontsize',fsz,...
        'xColor',grey,'yColor',grey,'color','k','ticklength',[0 0]);

% for all layers make contour plot
for iLay=1:gr.Nlay
    set(gca,'clim',hrange([1 end])); colorbar;
    
    xlabel(ax(iLay),'x [m]');
    ylabel(ax(iLay),'y [m]');
    title( ax(iLay), sprintf('head in layer %d',iLay))

    [c,hdl]=contourf(ax(iLay),gr.xm,gr.ym,H.values(:,:,iLay),hrange,'edgecolor','none'); % contour layer iLay
    clabel(c,hdl); % plot contour labels

    % plot wells in this layer
    Iw=find([well.iLay]==iLay);
    for iw=1:length(Iw)
        well(iw).plotXY('or','parent',ax(iLay));
    end

    % Plot drain
    if iLay==1, plot(ax(iLay),drn(:,1),drn(:,2),'color','g','linewidth',2); end
    
end

s={'Adapted example from MODFLOW 2000 manual, p89';''; 'Computed in mfLab';''; 'TO 2012-05-22'};
text(0.5,0.5,s,'HorizontalAlignment','center','color','yellow','fontsize',fsz);

%% Use showlayers

% h = showLayers(gr,'HK'); % ,well,'example');  % needs to be checked
%% Use zone budget to get budget overview

B=readBud([basename, '.bgt']); % read the Budget file

Zonearray=repmat(permute(1:gr.Nlay,[1,3,2]),[gr.Nlay,gr.Nx,1]);  % first create a zone array (here onze zone per layer)

zonebudget(B,Zonearray,1); % print budget of layer 1
zonebudget(B,Zonearray,2); % same for layer 2
zonebudget(B,Zonearray,3); % same for layer 3

zonebudget(B,Zonearray,[1 2 3]); % for all 3 layers totaled
