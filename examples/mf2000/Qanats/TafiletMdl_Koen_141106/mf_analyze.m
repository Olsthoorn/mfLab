%% Analyze Tafilalet model Coert Strikker, 2013_version Koen Geul
% Here the analyze of the model takes place. It has been updated to work
% with Koen_scenarios.m 
%
% KG 171014

fprintf('\n\n\n\n***** %s -- started *****\n\n',mfilename); 

%% Intro, loading of data

clear variables;   % start with clean workspace

load Scenarios     % load the scenarios
load name          % get name of model
load(basename)     % load model data save in mf_setup
load underneath    % get DEM etc, saved at end of mf_adapt

axProps = {'nextplot','add','xGrid','on','yGrid','on'}; % default figure porporties
formationColors = 'ymcg';
XLIM = gr.xGr([1 end]);
YLIM = gr.yGr([end 1]);

%% Read in computed head and cell by cell budget data

H = readDat([basename '.HDS']); % get the heads
B = readBud([basename '.BGT'],'t',[H.totim]); % read the budget file

%% Save the last head to be used as initial heads in a subsequent run
% SEE use in mf_adapt
if S.UseLH
    STRTHD   = H(end).values(:,:,1);   % #ok
    strthd   = gr.Z(:,:,end) + 0.5 * abs(diff(gr.Z(:,:,([1 end])),1,3));
    STRTHD(isnan(STRTHD)) = strthd(isnan(STRTHD));
    xc = gr.xc;
    yc = gr.yc;
    save INITIALHD STRTHD xc yc
end

%% Choose a cross section and see if there are khettaras intersecting this crosssection
% Choose your own cross section, and see whether khettaras are intesecting
% and in which layer they are intersecting

if S.KhetIntersect
    Khet_intersection(H,gr,khettaras);
end

%% Contour the heads
if S.ContourDEM
    figure('name','model overview','pos',screenPos(0.75),'renderer','zbuffer');

    axes(axProps{:},'xlim',XLIM,'ylim',YLIM);  % create axis, using props and limits
    xlabel('x UTM [m]'); ylabel('y UTM [m]');
    title('Heads in first layer at end of simulation period');

    % plot image
    A = imread('modelAreaImage.png'); image(gr.xGr,gr.yGr,A);

    % contour heads
    hrange= 750:1:950;
    [c,h]=gr.contourf(mean(H(end).values,3),hrange); clabel(c,h,'color','w');

    set(gca,'clim',hrange([1 end]));
    h = colorbar;  set(get(h,'title'),'string','Head [m wgs84]');
    set(gcf, 'renderer', 'zbuffer');
end

%% Isohypse map vs. situation 1959-1969
if S.Map
    figure('name','Isohypsen vs situation 2000-2010','pos',screenPos(0.75),'renderer','zbuffer'); hold on;
    A = imread('IsohypsesChap44p407.png');

    nIm = [31.582908 31.355736];
    eIm = [-4.504575 -4.296683]; 

    xIm = [357238 376668];
    yIm = [3495188 3469754];

    image(xIm,yIm,A);
    axis equal
    axis tight

    [c,h] = contour(gr.xm,gr.ym,H(end).values(:,:,1),-700:10:1000,'lineWidth',3);
    clabel(c,h);
end

%% HEAD relative to ground surface
if S.Head
    figure('name','head relative ground surface','pos',screenPos(0.75),'renderer','zbuffer');
    axes(axProps{:},'xlim',XLIM,'ylim',YLIM);
    xlabel('x UTM [m]'); ylabel('y UTM [m]');
    title('Head relative to ground surface. Pos and Neg separated with a thin white line');

    [c,h] = gr.contourf(H(end).values(:,:,1)-gr.Z(:,:,1)  ,-40:1:25);  clabel(c,h,'color','y');

    % Add white contour showing the zero line where head is at ground surface
    gr.contour (H(end).values(:,:,1)-gr.Z(:,:,1) ,[0 0],'w','lineWidth',1);

    hb = colorbar; set(get(hb,'title'),'string','\Phi - DEM [m]');
    set(gcf, 'renderer', 'zbuffer');
 
%      plot(gr.xm(89),gr.ym(24),'ko','markersize',10,'linewidth',3)

    %% DEPTH water in aquifer (head - bottom of aquifer)
    figure('name','water depth','pos',screenPos(0.75),'renderer','zbuffer');   
        % openGL instead of zBuffer necessary for transparency
    axes(axProps{:},'xlim',XLIM,'ylim',YLIM);  % create axis, using props and limits

    xlabel('x UTM [m]'); ylabel('y UTM [m]');
    title('Depth of groundwater above bedrock (=water depth)');
    
    hrange = 0.01:1:40;
    gr.contourf((H(end).values(:,:,1)-gr.Z(:,:,end)),hrange);
    set(gca,'clim',hrange([1 end]));

    hb = colorbar; set(get(hb,'title'),'string','water depth [m]');
    set(gcf, 'renderer', 'zbuffer');
    modelArea.plot
end

%% Plot transverse cross-sections
% Get the coordinates of the desired arbitrary profiles
if S.Cross
    profiles = lineObj( gr,[XLSData, 'JorfData.xls'],'profiles','NIL',gr.Z(:,:,1));

    % n profiles per page
    nfig = ceil(numel(profiles)/3);
    k=0; n=2;
    for ifig=1:nfig
        figure('name',sprintf('cross sections %d %d',k+1,k+n),...
               'pos',screenPos(0.75));
        xlim = [0 max([profiles(k+1:min(k+n,numel(profiles))).length])];

        for ipr=n:-1:1
            if ipr+k>numel(profiles), break; end
            subplot(n,1,ipr,axProps{:},'xlim',xlim);
            xlabel('s [m] along profile');
            ylabel('elevation [m]');
            title(['XSection along ' profiles(ipr+k).name]);
            profiles(ipr+k).fill2(gr,formationColors);
            profiles(ipr+k).plot2(gr,H(end).values,'linewidth',3); % heads
            if strcmpi('Guedima',profiles(ipr+k).name)
                profiles(ipr+k).plot2('k','linewidth',2);
            end
        end
        k=k+3;
    end
end

%% Draw profiles through the khettaras and also show the khettaras
if S.CrossKhet
    % we do 4 khettaras per page
    if exist('khettaras','var')
        Nk    = numel(khettaras);
        Npicture = 3;
        Npage = ceil(Nk/Npicture);
        xlim = [0 max([khettaras.length])];
        for ipage = 1:Npage
           figure('name',sprintf('khettaras page %d',ipage),'pos',screenPos(0.75));
           for j = 1:Npicture
               i = (ipage-1)*Npicture+j;
               if i>numel(khettaras), continue; end
               subplot(Npicture,1,j,axProps{:},'xlim',xlim);
               xlabel('s [m]');
               ylabel('elev [m]');
               title(khettaras(i).name);
               khettaras(i).fill2(gr,formationColors);
               khettaras(i).plot2('k','lineWidth',2);
               khettaras(i).plot2(gr,H(end).values,'b','linewidth',2);
           end
        end
    end
end

%% Draw pofiles of all lines and areas that are in your models
if exist('headBoundaries','var') && exist('fluxBoundaries','var') && exist('riversDRN','var')
    boundaries = [headBoundaries(:); fluxBoundaries(:); riversDRN(:)];
    Nb = numel(boundaries);
    Npicture = 2;
    Npage = ceil(Nb/Npicture);
    xlim = [0 max([boundaries.length])];
    for ip = 1:Npage
        figure('name',sprintf('boundary profiles page %d',ip),'pos',screenPos(0.75));
        for i = 1:Npicture
            ib = (ip-1)*Npicture+i;
            if ib>Nb, continue; end
            subplot(Npicture,1,i,axProps{:},'xlim',xlim);
            xlabel('s [m]'); ylabel('elevation [m]');
            title(sprintf('%s',boundaries(ib).name));
            boundaries(ib).fill2(gr,formationColors);
            boundaries(ib).plot2('k');
            boundaries(ib).plot2(gr,H(end).values,'b','linewidth',2);
        end
    end
end

%% PROFILES AROUND pumpAreas
if S.ProfilePump
    % This given a profile like running around the area along its circumference
    % Looks strange, is not very useful but is correct.
    if exist('pumpAreas','var') % #ok
        figure('name','pumpAreas profiles along their circumference','pos',screenPos(0.75));
        xlim = [0 max([pumpAreas.length])];
        Na = numel(pumpAreas);
        for ia = 1:Na
            subplot(Na,1,ia,axProps{:},'xlim',xlim);
            xlabel('s [m]'); ylabel('elevation [m]');
            title(sprintf('%s',pumpAreas(ia).name));
            pumpAreas(ia).fill2(gr,formationColors);
            pumpAreas(ia).plot2('k');
            pumpAreas(ia).plot2(gr,H(end).values,'b','linewidth',2);
        end
    end
end

%% GEOLOGIC XSection in 3D with labels
if S.Cross
    demContour     = lineObj( gr,[XLSData 'JorfData.xls'],'demContour','NIL',gr.Z(:,:,1));

    Npr  = numel(profiles); % number of profiles
    figure('name','geological cross sections in 3D','pos',screenPos(0.75)); % a new figure

    axes(axProps{:},'xlim',XLIM,'ylim',YLIM,'zgrid','on');
    xlabel('x [m]'); ylabel('y [m]');

    names = {profiles.name};
    for j= ceil(numel(names)/5):-1:1
        J = 5*(j-1) + (1:5); J=J(J<numel(names));
        ttl{1,j} = sprintf(' %s',profiles(J).name);
    end
    title(['profiles:',ttl]);

    for ipr=Npr:-1:1 
        profiles(ipr).fill3(gr,formationColors);             % profile
        profiles(ipr).plot3(gr,H(end).values,'lineWidth',2); % heads
        profiles(ipr).label3;
    end

    modelArea.plot3('r','lineWidth',1);
    if exist('outcrops','var')
        outcrops.plot3('k','lineWidth',1);
    end

    demContour.plot3('k','lineWidth',1);

    view(60,45)
    view(2)
end

%% SHOW WHAT IS IN MY MODFLOW STRESS FILES
if 0
    figure('name','check head boundaries','pos',screenPos(0.75)); hold on;
    xlabel('x UTM [m]'); ylabel('y UTM [m]'); title('fixed-head boundaries');
    [HB] = headBoundaries.plot('ro','lineWidth',2);

    d = dir([basename '.*']);
    for i=numel(d):-1:1
        [~,~,ext] = fileparts(d(i).name);
        if ~ismember(ext(2:end),{'WEL','DRN','GHB','RIV','CHD'})
            d(i)=[];
        end
    end
    leg{numel(d)}='';
    for i=1:numel(d)
        leg{i} = d(i).name;
        gr.showWhatIsInMyStressFile(d(i).name,[mf_color(i+1) 'o']);
    end
    legend(leg{:});
end

%% PLOt HEADS ON A 3D SURFACE
% nice but not necessary very usefull

if S.DView
    figure('pos',screenPos(0.75)); % #ok
    hold on; grid on;
    xlabel('x [m]'); ylabel('y [m]');
    title('Heads in the cells as a 3D surface');
    surf(gr.xGr,gr.yGr,gr.ZGR(:,:,1),H(end).values(:,:,1));
    view(3);
end

%% PLOT TIME SERIES FOR HEADS IN PIEZOMETERS (computed)
% If the piezom objects are loaded with their measured heads. They too
% could be plotten, together with the computed ones as is done here.

if S.Piezom
    if 0
         hd1        = piezom.plotHead(gr,H,false,'lineWidth',2); %#ok
        [hd2,names] = piezom.plotHead(gr,H,true ,'lineWidth',2);
    else
        piezom           = piezom.plotHead(gr,H,false,'meas','lineWidth',2);
        [PERnams,PERvals]= getPeriods(basename,'PER');

        P1028 = PERvals(:,strmatchi('P1028',PERnams)); 
        hTime = [H.time];
        hTime = hTime(P1028>100);
        P1028 = P1028(P1028>100);
        plot(hTime,P1028,'ro');   

        P1029 = PERvals(:,strmatchi('P1029',PERnams)); 
        hTime  = [H.time];
        hTime    = hTime(P1029>100);
        P1029 = P1029(P1029>100);
        plot(hTime,P1029,'gx','markersize',10);   

        P3628 = PERvals(:,strmatchi('P3628',PERnams)); 
        hTime  = [H.time];
        hTime    = hTime(P3628>100);
        P3628 = P3628(P3628>100);
        plot(hTime,P3628,'m+');   

        piezom = piezom.plotHead(gr,H,true ,'meas','lineWidth',2);
    end
end

%% DISCHARGES of Khettaras
if S.QKhettaras
    dimension = 'l/s';  % also possible 'm3/d', 'm3/h', 'm3/y' or 'm3/a'  'Mm3/y' or 'Mm3/a'
    printAllObjDischarges;
end

%% Flow across section
if S.Cross
    % BE AWARE: THIS TAKES A LONG TIME if done for many long cross secions and many time steps!!!!

    % Therefore only a few interesting cross sections are computed
    % and this for only a few stress periods, (i.e. every 24 stress periods,
    % which implies every 2 years (see time) above list)

    % You can choose dimensions 'Mm3/a m3/d l/s m3/h ...

    % printQacrossProfiles;  % this does everything, don't do it because it may
                             % take half an hour perhaps

    % This is better, choose your cross sections and your stress periods and also a suitable dimension:
    profiles(strmatchi({'AA-Fezna outcrop','Fezna outcrop - cretaceous', 'outflow AA-Krayr', 'Cross04','Catchments'},{profiles.name})).Qacross(gr,B(24:24:end),'l/s');

    % with spec, the value is divided by the cross section length (specific
    % discharge instead of total discharge. Therefore also teh dimension has
    % been chosen differently to get interpretable numbers.
    profiles(strmatchi({'AA-Fezna outcrop','Fezna outcrop - cretaceous', 'outflow AA-Krayr', 'Cross04','Catchments'},{profiles.name})).Qacross(gr,B(24:24:end),'m3/d','spec');

    % profiles(strmatchi({'Cross01', 'Cross04', 'Cross10', 'Cross13'},{profiles.name})).UserData
end

%% Compute the flux for all boundaries

% This is the flux from the budget file for the boundary's type (WEL, DRN
% etc) at total and per cell. See the output struct and help lineObj/Qtype.

% qF = fluxBoundaries.Qtype(B); % total flow by fluxBoundaries for all objects and stress periods
% qH = headBoundaries.Qtype(B);

%% Compute the total transmissivity of a profile
% The transmissivity is compute as the integral of the sum of the wet thickness of the
% layers times their horizontal conductivity and lengtn

% kD = integral( sum_over all layers(kh * Dwet) ds (where s is taken along the length of the profile
% The dimension is [L3/T] and not [L2/T], it is in fact a conductance because of the integration
% along its length.

if S.Cross
    [kD,profiles] = profiles.kD(gr,H,HK); % Notice, kD can be any lineObj or area2Obj
end

%% CHECKING RESULTS:
%% Do the heads along fixed-head boundaries coincides with the input.
%% Are the flows consistent with the input for all boundaries with given flows.
%% Does the total recharge match the input?
%% Calculate by hand (in Matlab) the total inflow and outflow form specific boundaries like drains, riversDRN and general head boundaries.
%% Does the total storage match with the change of heads?
%% Does the model results conform to the modeler?s expectations? Especially, is the discharge from the Khettaras in accordance with field measurements?
%% Does the computed head agree with known heads in piezometers?

fprintf('***** %s -- finished. Running it took %g seconds *****\n',mfilename,toc);


    