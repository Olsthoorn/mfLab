%% Analyze Tafilalet model Coert Strikker, 2013

%% Intro, loading of data

close all;         % close all previous pictures
clear variables;   % start with clean workspace

load name          % get name of model
load(basename)     % load model data save in mf_setup
load underneath    % get DEM etc, saved at end of mf_adapt

% some colors
Blueish   =[ 18/255,104/255,179/255];
Reddish   =[237/255, 36/255, 38/255];
Greenish  =[155/255,190/255, 61/255];
Purplish  =[123/255, 45/255,116/255];
Yellowish =[  1    ,199/255,   0   ];
LightBlue =[77/255 ,190/255,238/255];



xlim = gr.xGr([1 end]); ylim= gr.yGr([end 1]); % figure limits
axProps = {'nextplot','add' };                  % default figure porporties

H = readDat([basename '.HDS']);               % get the heads
B = readBud([basename '.BGT'],'t',[H.totim]); % read the budget file, add time

%% Contour heads

figure('pos',screenPos(0.75)); % create figure to plot on

iax = 1;

ax(iax) = axes(axProps{:},'xlim',xlim,'ylim',ylim);  % create axis, using props and limits

xlabel(ax(iax),'x [m]');
ylabel(ax(iax),'y [m]');
title( ax(iax),'Heads in first layer at end of simulation period');

hrange= ContourRange(H,100);   % Get a bunch of contour elevations.

contourf(ax(iax),gr.xm,gr.ym,mean(H(end).values,3),hrange); % contour in full color

gr.plotGrid(ax(iax),'c');  % cyan (c) is default grid color

hb= colorbar; set(get(hb,'title'),'string','head');

%% Plot objects and label them group by group
% It can also be done in a general loop as was shown in the Jorf case,

modelArea.plot( ax(iax),'color',[0.3 0.3 0.3])    ; modelArea.label(ax(iax));
outcrops.fill(ax(iax),'m')              ; outcrops.label( ax(iax));
pumpAreas.plot(ax(iax),'g')             ; pumpAreas.label(ax(iax));
rivers.plot(ax(iax),'b','lineWidth',3)  ; rivers.label(ax(iax));
khettaras.plot(ax(iax)),'k'             ; khettaras.label(ax(iax));
headBoundaries.plot(ax(iax),'color',Blueish); headBoundaries.label(ax(iax));
fluxBoundaries.plot(ax(iax),'color',Greenish ); fluxBoundaries.label(ax(iax));
genHeadBoundary.plot(ax(iax),'gp-.')    ; genHeadBoundary.label(ax(iax));
wells.plot(ax(iax),8,'k')               ; wells.label(ax(iax));
piezom.fill(ax(iax),8,DEM)              ; piezom.label(ax(iax));
profiles.plot(ax(iax),'color',Greenish) ; profiles.label(ax(iax));

%% Show the profiles in flat sections (2D projected)

Np = numel(profiles); % number of profiles
Na = numel(ax);       % number of currently active axes
figure('pos',screenPos(0.75)); % a new figure

for ipr=Np:-1:1   % running backwards, so no preallocation is needed
                  % this informs Matlab how many upfront
    iax = ipr+Na;   % axis number
    ax(iax) = subplot(Np,1,ipr,'nextPlot','add','xlim',[0 1200],'xgrid','on','ygrid','on');
    xlabel(ax(iax),'s [m] along profile');
    ylabel(ax(iax),'elevation [m]');
    title( ax(iax),['XSection along ' profiles(ipr).name]);
    profiles(ipr).fill2(ax(iax),gr.Z,'ymc');          % profile
    profiles(ipr).plot2(ax(iax),H(end).values); % heads
end

%% Do the khettaras have a constant inclination --> project them on a vertical plane through them
figure('pos',screenPos(0.75));
iax=length(ax)+1;
ax(iax) = axes('nextPlot','add'); xlabel('s along section [m]'); ylabel('elevation [m]'); title('khettaras');
khettaras.fill2(ax(iax),gr.Z,'ymc');
khettaras.plot2(ax(iax));


%% Show the profiles in 3D

Np = numel(profiles); % number of profiles
Na = numel(ax);       % number of currently active axes
figure('pos',screenPos(0.75)); % a new figure

iax = 1+Na;   % axis number
ax(iax) = axes('nextPlot','add','xlim',xlim,'ylim',ylim,'xgrid','on','ygrid','on','zgrid','on');

xlabel(ax(iax),'x UTM [m]');
ylabel(ax(iax),'y UTM [m]');

title( ax(iax),['Profiles in 3D space: ', sprintf(' %s',profiles.name)]);

for ipr=Np:-1:1   % running backwards, so no preallocation is needed
                  % this informs Matlab how many upfront
    profiles(ipr).fill3(ax(iax),gr.Z,'ymc');          % profile
    profiles(ipr).plot3(ax(iax),H(end).values); % heads
end

modelArea.plot3(ax(iax),'r','lineWidth',3);
outcrops.plot3( ax(iax),'k','lineWidth',3);

%% Plot the heads on a surface in 3D
% Each cell will be placed in 3D space (i.e. on a surface) at its proper
% elevation and colored according to its head. All cells outside the model
% area, i.e. where IBOUND=0 will not be colored.

figure('pos',screenPos(0.75)); hold on; grid on;
xlabel('x [m]'); ylabel('y [m]');

iax     = numel(ax)+1;
ax(iax) = gca;

title(ax(iax),'Heads in the cells as a 3D surface');

surf(gr.xGr,gr.yGr,gr.ZGR(:,:,1),H(end).values(:,:,1),'parent',ax(iax));
khettaras.plot3('k');

hb = colorbar; set(get(hb,'title'),'string','elevation');

view(3);

%% Getting the discharges

% fprintf('Discharge of all khettaras in a list:\n');
% khettaras.printQ(B);
%
% To get the cell discharges 
% if ~isempty(khettaras(1).UserData) && isfield(khettaras(1).UserData,'Qcell')
%     khettaras(1).UserData.Qcell
%     fprintf('\n');
% end
%
% To get the specific discharges
% if ~isempty(khettaras(1).UserData) && isfield(khettaras(1).userData,'Qspec')
%     khettaras(1).UserData.Qspec  
%     fprintf('\n');
% end

%% Water balances using the stress objects

%% get the flows and fluxes of all stresses
fprintf('Q of headBoundaries\n');   headBoundaries = headBoundaries.printQ(B);
fprintf('Q of pumpAreas\n');        pumpAreas      = pumpAreas.printQ(B);
fprintf('Q of khettaras\n');        khettaras      = khettaras.printQ(B);
fprintf('Q of rivers\n');           rivers         = rivers.printQ(B);
fprintf('Q of genHeadBoundary\n');  genHeadBoundary = genHeadBoundary.printQ(B);
fprintf('Q of wells\n');            wells          = wells.printQ(B);
fprintf('Q of fluxBoundaries\n');   fluxBoundaries = fluxBoundaries.printQ(B);

%% Total water balance of all stress objects combined
% the sum adds them line by line.
% vertcat is required if there are more then one object, e.d. there area 8
% wells. sum(vertcat(wells.Q)) adds their flows, then the overall sum adds
% them with those of the other objects.

S = sum(...
     [headBoundaries.Q;...
      rivers.Q;...
      genHeadBoundary.Q;...
      sum(vertcat(wells.Q));
      sum(vertcat(khettaras.Q));...
      sum(vertcat(fluxBoundaries.Q));...
      pumpAreas.Q])

fprintf('\n\nThe overall water budget of all the stresses in the list reads: <<<<<%g m3/d>>>>>\n',sum(S));

  
%% Done