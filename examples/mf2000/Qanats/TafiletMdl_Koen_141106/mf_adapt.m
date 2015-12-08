%% Construction of Tafilalt model Coert Stikker Juni 2013_Version Koen Geul
% Here the all the values for the model are obtained. It has been updated
% to work with Koen_scenarios.m.
%
% KG 171014

load Scenarios % scenarios form koen_scenarios

fprintf('\n\n\n\n***** %s -- started *****\n\n',mfilename); tic;

%% NAME and DATA LOCATIONS of this project
% try to find a previous waterbalance, if not, it will not be used
try
    load Waterbalance
    fprintf('.. Waterbalance found!\n');
catch
    fprintf('.. Waterbalance not found (this can happen in the first run)\n');
end

basename= 'Jorf';
save ('name','basename');

XLSData = ['.' filesep 'XLSData'   filesep]; % EXCEL files directory
MATFiles = ['.' filesep 'MATFiles' filesep]; % Storage files directory
MFiles  = ['.' filesep 'MFiles' filesep];      % Storage files directory

addpath(MATFiles); % Add to path so they can be located
addpath(MFiles);  % Add to path so they can be located

esriLayers = {'DEM','silt','gravel','limestone','conglomerate'}; 
    % names of layer in data files

MINDZ = 0.01; % min allowed layer thickness

%% HYDRAULIC LAYER PARAMETERS

hk   = [1,260,40,86]; % [m/d] hor conductivity of the 3 geological ESRI layers
S.sy = 0.25;          % [ - ] specific yield for transient simulations
ss   = 1e-5;          % [1/m] elastic storativity

%% Loading existing grid file, else create new grid

if exist('gr.mat','file')
   load gr
else
    %% GET DEM
    % Here we get the dem directly from strm_36_06_Jorf_???m.mat in folder ./DEM_SRTM
    % Generating new DEMs can be done with the mfile in that directory.
    % Choose from 4 resolution by changing the value before the m in oneof 090 180 270 or 360
    load(fullfile('DEM_SRTM','strm_36_06_Jorf_270m.mat'));

    %% Generate the model GRID using the coordinates that come along with the loaded DEM

    % We don't know the layer elevations yet, so use use [0 -1] instead as dummy
    % We will change this after we have the layer thicknesses.

    gr = gridObj(xGr,yGr,[0 -1],'MINDZ',MINDZ);

    % Layer thickness by interpolating the Google Earth digitized contours

    Z = zeros(gr.Ny,gr.Nx,numel(esriLayers));
    Z(:,:,1) = DEM;
    for iLay=2:numel(esriLayers)  % notice that we start at 2 (skip DEM).
        % The interpolatin is done in this function
        D = gr.interpFromKMLpaths(fullfile('kmlData','geology',esriLayers{iLay}));
        D(isnan(D) | D<MINDZ)= MINDZ;  % remove NaN's and assert MINDZ 
        Z(:,:,iLay) = Z(:,:,iLay-1)-D; % elevation of current layer bottom
    end

    % Here comes our final gridObj with all layer elevations included.
    gr = gridObj(gr.xGr,gr.yGr,Z,'MINDZ',gr.MINDZ);

    save gr gr
end

%% BOUNDARIES OF THE DIFFERENT REGIONS
% Get boundary of model and of outcrops, so that we can cutout the cells
% that do not take part in the mdoel.

modelArea = area2Obj(gr,[XLSData,'JorfData.xls'],'extent','NIL','name','model');
outcrops  = area2Obj(gr,[XLSData,'JorfData.xls'],'extent','NIL','name','out');

% Cutting out the unwanted model cells is done with IBOUND in full 3D, but
% here we do this in 2D and apply it later also to the other deeper layers.
ibound= zeros(gr.Ny,gr.Nx);               % top layer of IBOUND
ibound(vertcat(modelArea.Idx)) = true;    % active model area

if exist('outcrops','var');
    % cutout each outcrop, by setting the corresponding cells to false in
    % ibound (false and 0 are treated the same here).
    for ioc=1:numel(outcrops), ibound([outcrops(ioc).Idx]) = false; end
end

%% Some small mistakes in ibound
% there are some small mistakes in the ibound boundaries. This is usefull
% for EVT. Otherwise EVT will not work. At these point the bottom of the
% aquifer is too close to the surface ( < S.DSurf). These points need to be
% excluded, otherwise EVT does not work.

SmallLay = find( sum(gr.DZ,3) < (S.DSurf + 0.50) & ibound > 0 );
ibound(SmallLay) = false;

%% IBOUND
% Exclude any area outside the model and within the outcrops. We already
% created the first layer of the ibound. Here it will be applied to all the
% layers.

IBOUND = bsxfun(@times,ibound,ones(1,1,gr.Nlay));

%% Contour the DEM + Layers
if S.LayerDEM
    figure('name','DEM','position',screenPos(0.75)); hold on;
    title('Tafilalt Model, Coert Strikker, 2013, DEM (SRTM 3'''' about 90 m)');
    xlabel('utm x [m]'); ylabel('utm y [m]');

    gr.contourf(gr.Z(:,:,1) .* ibound, 750:2:1200); colorbar;
    axis('equal'); axis('tight');
    hb = colorbar; set(get(hb,'title'),'string','elevation');
    set(gcf, 'renderer', 'zbuffer');

    % THICKNESS OF EACH LAYER
    for iLay=1:gr.Nlay
        figure('name',sprintf('D of layer %d',iLay),'pos',screenPos(0.75));
        axes('nextplot','add');
        title(sprintf('Thickness of Layer %s',esriLayers{iLay+1}));
        xlabel('x [m]'); ylabel('y [m]');
        [c,h] = gr.contourf(gr.DZ(:,:,iLay) .* ibound, 0:1:50); clabel(c,h,'color','w');
        h = colorbar; set(get(h,'title'),'string','Layer thickenss [m]');

        % Additionally plot the Google Eartch digitized contour lines themselves
        % on the obtained thicknes map 
        Dlines = dir(fullfile('kmlData','geology',[esriLayers{iLay+1} '*.kml']));
        if isempty(Dlines)
            error(['Can''t find files <<%s>>\n',...
                   'Remdy: 1) Check the spelling including capital and lowercase !\n',...
                   '       2) Check that you run mf_adapt from the correct directory.'],...
                   fullfile('kmlData','geology',[formationName '*.kml']));
        end

        for j=1:numel(Dlines)
            [xC,yC] = kmlPath2UTM(fullfile('kmlData','geology',Dlines(j).name));
            plot(xC,yC,[mf_color(j) '--'],'lineWidth',2);
        end


        axis('equal'); axis('tight');
        hb = colorbar; set(get(hb,'title'),'string',sprintf('D of %s',esriLayers{iLay+1}));

    end
end

%% THE RECHARGE, probably the most important part in the model
% We compute the recharge based on the monthly precipitation by
% considering the monthly precipitation as a single shower, subtrating its
% initial loss and multiplying the remainder with a recharge coefficient.
% These properties are specific for each catchment. Downstream catchments
% add runoff from upstream cathcmentts to their precipitation before
% subtracting the threshold.

    % Get the catchments directly form the kmlFolderFile
    catchment        = kmlPathsObj(fullfile('kmlData','areas','catchments'));
    catchment(end+1) = kmlPathsObj(fullfile('kmlData','areas','modelcontour'));
    catchment(end).name = 'restZone';

    NCatch = numel(catchment);

    % "restment" is the model area that is not part of any of the defined
    % catchments

    % Get catchment properties from the worksheet 'catchments'
    % User: Make sure the list in the worsheet matches with the names of 
    % the catchments contained in the kml file
    [CatchmHdr,CatchVals,CatchmTxtHdr,CatchTxt] = ...
        getExcelData([XLSData 'jorfData'],'catchments','hor');

    % Extract the catchment properties
    initialLoss  = CatchVals(:,strmatchi('initialLoss',CatchmHdr));    
                   % inital loss of this cathcment [m]
    rf           = CatchVals(:,strmatchi('rf'         ,CatchmHdr));    
                   % runoff factor for this catchment 
    gf           = CatchVals(:,strmatchi('gf'         ,CatchmHdr));    
                   % Ground flow for this catchment 
    S0           = CatchVals(:,strmatchi('S0'         ,CatchmHdr));    
                   % Initial Storage for this catchment 
    name         = CatchTxt (:,strmatchi('name'       ,CatchmTxtHdr)); 
                   % catchment name
    downStr      = CatchTxt (:,strmatchi('downStr'    ,CatchmTxtHdr)); 
                   % downstream catchment name

    % Add catchment properties to catchments, so that each knows its own
    % properties.
    for ic=1:NCatch
        Iupstr = strmatchi(name(ic),downStr);
        if ~Iupstr
            catchment(ic).UserData.upstr = {};
        else
            catchment(ic).UserData.upstr = name(Iupstr);
        end

        % if recharge is selected it will assign values, oherwise it will be zero
        % some catchments are selected not to contribute to the recharge by
        % a high initial loss, those values are higher than 0.05. by if you
        % can select them and assign the orginal value. for the other
        % values you have S.InitialLoss, which is defined in scenarios.
        % Same count for the runoff factor, some have a different one, this
        % value will remain the same.
        
        if initialLoss(ic) < 0.05
            catchment(ic).UserData.initialLoss = S.InitialLoss;
        else
            catchment(ic).UserData.initialLoss = initialLoss(ic);
        end
        
        if rf(ic) > 0.7
            catchment(ic).UserData.rf = S.RunOff;
        else
            catchment(ic).UserData.rf = rf(ic);
        end
        
        catchment(ic).UserData.gf          = gf(ic);
        catchment(ic).UserData.S0          = S0(ic);
        catchment(ic).UserData.visited     = false;
    end

    if S.SpyCatch 
        % used for control, only if S.Recharge = true.
        for ic=1:numel(catchment)
        fprintf('%15s receives water from <<%s >>\n',...
            catchment(ic).name,...
            sprintf(' %12s',catchment(ic).UserData.upstr{:}));
        end
    end

    % Get the PRECIP (not RECH) from PER sheet and translate it into recharge
    % catchments which receive also water from upstream catchments (monthly values !!)
    [PERnams,PERvals,NPER] = getPeriods(basename,'PER',{'PERLEN','Precip'});
    if numel(S.PercentRecharge)==1
        RechFac = S.PercentRecharge;  % 12 months times 11 years
    else
        RechFac = (S.PercentRecharge * ones(1,12)).';  % 12 months times 11 years
    end
    PERvals(:,1) = PERvals(:,1) .* RechFac(:);

    precip = PERvals(:,strmatchi('precip',PERnams));  % precip per stress period (month)
    perlen = PERvals(:,strmatchi('PERLEN',PERnams));  % length of stress period (days)

    % Recursively get the runoff and infiltration for all catchments and runoff
    % areas (also treated as catchments). The runoff and rcharge for all areas
    % is stored in the catchments's UserData as fields.
    catchment = getRech(catchment,precip,perlen,S,'balance','all','Lines',10);

    % Generate the recharge array, with recharge for each model cell and each
    % stress period
    RECH      = zeros(gr.Ny,gr.Nx,NPER);

    % We already computed the recharge and runoff for all catchements above in
    % the function getRech(catchments). Now distribute this recharge over the
    % cells of the model with each catchment.

    % RestZone is a zone array with true for every cell within the model
    % contour and false outside the contour. We cutout the cells that are
    % within the the model but not in one of the catchments. This way each
    % of the catchments including the restzone will obtain its own recharge.
    catchment(NCatch).UserData.In = inpolygon...
        (gr.Xm,gr.Ym,catchment(NCatch).X,catchment(NCatch).Y);
    
    % The model cells pertaining to the catchments will be shown using spy()
    if S.SpyCatch
        figure('name','catchments');
        xlabel('pixels'); ylabel('pixels');
        title('catchments');
        spy(catchment(NCatch).UserData.In); hold on;
    end
    
    % Make sure the In field of the catchments corresponds excludes restZone
    for ic=NCatch:-1:1
        % Cells in this catchment:
        inZone=inpolygon(gr.Xm,gr.Ym,catchment(ic).X,catchment(ic).Y);
        catchment(ic).UserData.In = inZone;

        if S.SpyCatch
            spy(inZone,mf_color(ic,'rgmcky'));   % inspect if correct
            set(gcf, 'renderer', 'zbuffer');
        end

        % Add the cathment recharge to the cells of the model within the respective
        % catchments (these valuesa are always in m/month)

        if ic ~= NCatch
            % cutout this catchment from restZone    
            catchment(NCatch).UserData.In = ...
                catchment(NCatch).UserData.In & ~inZone;
        end
    end

if S.Recharge
    % Finally fill RECH with the proper zone values
    for ic=1:numel(catchment)
        RECH = RECH + bsxfun(@times,catchment(ic).UserData.In,...
            XS(catchment(ic).UserData.Rch));
    end
else
    % if S.Recharge = false; RECH will only be 0.
    RECH = zeros(size(gr.AREA3));
end

%% Check how much recharge we got
if S.Recharge
    fprintf('Yearly recharge per catchment in mm:\n');
    for ic = 1:numel(catchment)
        fprintf('%20s %8.1f mm/y',catchment(ic).name,1000 * mean(catchment(ic).UserData.Rch) * aYear);
        if ~isempty(catchment(ic).UserData.upstr)
            fprintf('. It receives runoff from {%s }.',sprintfs(' %s',catchment(ic).UserData.upstr)); 
        end
        fprintf('\n');
    end

    fprintf('%20s %8.1f mm/y\n','Notice: precip = ',1000 * mean(precip(:)./perlen) * aYear);
    fprintf('Small differences may occur if catchments have overlaps due to inaccurate digitization.\n');
    fprintf('This may or may not be acceptable.\nI think it is acceptable if small compared to the total recharge.\n');
end


%% EVT
% Here the parameters for the EVT are created. If ~S.EVT it will all be 0.
if S.EVT 
    [PERnams,PERvals,NPER] = getPeriods(basename,'PER',{'EVTR'});
    PERvals(:,1) = PERvals(:,1) * S.PercentEVT; 
       
    evtr   = PERvals(:,strmatchi('EVTR',PERnams));    % Evaporation per stress period
    
    SURF = max( gr.Z(:,:,1) - S.DSurf,gr.Z(:,:,end));
    EXDP = SURF - max( SURF - S.ExDepth,gr.Z(:,:,end));
    
    SURF = bsxfun(@times,SURF,ones(1,1,NPER));
    EXDP = bsxfun(@times,EXDP,ones(1,1,NPER));
    
    Aevap = ones(size(gr.AREA)); % create top layer with all zeros
    
    EVTR = bsxfun(@times,Aevap,XS(evtr)); % calculate EVTR per cell
else
    SURF = zeros(gr.Ny,gr.Nx,NPER);
    EXDP = zeros(gr.Ny,gr.Nx,NPER);
    EVTR = zeros(gr.Ny,gr.Nx,NPER);
end

%% GENERATE and ADAPT the HK (horizontal k values)
HK    = gr.const(hk);

% Silt factor is used the increase the conductivity of the silt layer
% towards the mountain faces as generally it gets higher the steeper the
% terrain.
siltFactor = gr.interpFromKMLpaths(fullfile('kmlData','geology','hkSilt'));
HK(:,:,1) = HK(:,:,1).*siltFactor;

%% IF WE WANT TO COMBINE LAYERS IT MUST BE DONE HERE!!

% For convergence reasons and faster testing of the model, it is useful to
% combinine layers as very thin layers often cause issues. This is done
% here using the simpel fromTO array below, telling which layers of the
% original model will be included in which layers from the new model.
combineSomeLayers = true;

fromTo = [1 1;  % from is layers in orginal model, to in new model
          2 2;
          3 3;
          4 4];
      
% Here we carry out the coversion:
NlayFr = numel(unique(fromTo(:,1)));
NlayTo = numel(unique(fromTo(:,2)));
DZto   = zeros(gr.Ny, gr.Nx, NlayTo);
kDto   = zeros(gr.Ny, gr.Nx, NlayTo);

for iLay= 1:NlayFr % Fr s the from model To is the to mdoel
    iFr = fromTo(iLay,1);
    iTo = fromTo(iLay,2);
    DZto(:,:,iTo) = DZto(:,:,iTo)+gr.DZ(:,:,iFr);
    kDto(:,:,iTo) = kDto(:,:,iTo)+gr.DZ(:,:,iFr).*HK( :,:,iFr);
end

HK = kDto./DZto;   % the new HK
Zto= bsxfun(@times,gr.Z(:,:,1),ones(1,1,NlayTo+1)); % the new Z array

for iLay=1:NlayTo
    Zto(:,:,iLay+1) = Zto(:,:,iLay)-DZto(:,:,iLay);
end

gr = gridObj(gr.xGr,gr.yGr,Zto,'MINDZ',MINDZ); % the new grid

%% The other new model arrays (a bit simplistic, but ok here)
VK = HK/5;
SY = gr.const(S.sy);
SS = gr.const(ss);

%% STARTHD initial heads

% load a previously saved head layer as initial heads
msgId = 'mf_adapt:missingSTRTHD';
try
    % try to load file INITIALHD.mat, which, if present contains STRTHD (previously saved).
    load INITIALHD;
    if any(size(STRTHD(:,:,1))~=[gr.Ny,gr.Nx])        
        % if the size does not match that of the grid, because you have changed
        % the resolution of the model, then simply interpolate the old start
        % heads to the cell centers of the grid. xc and yc are contained in
        % INITIALHD.mat
        STRTHD = interp2(xc,yc,STRTHD(:,:,1),gr.xm,gr.ym);
    end

catch ME
    % if file does not exist, warn user, but continue (we will use a default)
    warning('on',msgId);
    warning(msgId,...
            ME.message,...
            ['No previously saved STRTHD found (looking for file <<INITIIALHD.mat>>\n',...
             'Will use default initial heads']);
     warning('off',msgId);
     
    % use default STRTHD
    STRTHD= bsxfun(@times,gr.Z(:,:,1)-0.5*sum(gr.DZ,3),ones(1,1,gr.Nlay));
end

% when you have a layer of start heads, apply it to all model layers
STRTHD = bsxfun(@times,STRTHD(:,:,1),ones(1,1,gr.Nlay));

%% OBJECTS (HYDROLOGIC FEATURES)

% Add any pointObj, lineObj and area2Obj of any stress type to the model

if S.ConstantH
    headBoundaries = lineObj( gr,[XLSData,'JorfData.xls'],'boundaries',...
        'CHD','type','head');
end

if S.GradientBS
    gradient = sqrt(median(median(diff(gr.Z(:,:,1),1,2)./diff(...
        gr.XM(:,:,1),1,2)))^2 + median(median(diff(gr.Z(:,:,1),1,1)./...
        diff(gr.YM(:,:,1),1,1)))^2);
    North = ismember({headBoundaries.name},'north');
    headBoundaries = [headBoundaries(North) ...
        headBoundaries(~North).setGradient(gr,gradient,HK)];
end

if S.River 
    % if S.River, it will assign values with the percentage, else it will 
    % all be 0
    riversQ   = lineObj( gr,[XLSData,'JorfData.xls'],'riversQ','RIV',...
        'type','riv');
    
    for iR = 1:length(riversQ)
        C = strmatchi('c',riversQ(iR).vertexHdr);
        facRiv = (S.PercentRiver * ones(1,12))';
        facRiv = mean(facRiv(:)); % recharge factor for the time being not yet splito over the years and the months
        riversQ(iR).vertex(:,C)       = riversQ(iR).vertex(1,C) .* facRiv(:);
        riversQ(iR).cellLineVals(:,C) = riversQ(iR).cellLineVals(1,C) .* facRiv(:);
    end
else
    riversQ = lineObj(gr,[XLSData,'JorfData.xls'],'riversQ','RIV',...
        'type','riv');
    
    % if c = 0, no water will infiltrate
    for iR = 1:length(riversQ)
        C = strmatchi('c',riversQ(iR).vertexHdr);
        riversQ(iR).vertex(:,C) = 0;
        riversQ(iR).cellLineVals(:,C) = 0;
    end
end

if S.Khettaras 
    % the outflow of the khettaras will be used for the irrigation of
    % Khettara water in the Palmaraies in the next loop
    khettaras = lineObj( gr,[XLSData, 'JorfData.xls'],'khettaras','DRN');
    pumpAreas = area2Obj(gr,[XLSData, 'JorfData.xls'],'Pumping','WEL');
    FirstI    = find(IBOUND,1); % Used to change Idx of khettaras to get 
                                % better plots. they will be located at an
                                % location where they will have no effect
                                % but cant be seen on the plots

    if ~S.SuperKhetFJ && ~S.SuperKhetHa && ~S.SuperKhetFezna
        Super = strmatchi('Super',{khettaras.name});
        for iK=1:length(Super)
            HKe = strmatchi('H',khettaras(Super(iK)).vertexHdr);
            Z = strmatchi('Z',khettaras(Super(iK)).vertexHdr);
            C = strmatchi('c',khettaras(Super(iK)).vertexHdr);
            

            % high values, so no water can flow into the drain (this is the
            % elevation),
            khettaras(Super(iK)).vertex(:,HKe)        = 10^10;    
            khettaras(Super(iK)).vertex(:,Z(2))       = 10^10;
            khettaras(Super(iK)).cellLineVals(:,Z(2)) = 10^10;
            khettaras(Super(iK)).cellLineVals(:,HKe)  = 10^10;
            khettaras(Super(iK)).vertex(:,Z(1))       = 10^10;
            khettaras(Super(iK)).cellLineVals(:,Z(1)) = 10^10;
            khettaras(Super(iK)).vertex(:,C(1))       = 10^10; 
            khettaras(Super(iK)).cellLineVals(:,C(1)) = 10^10;
            khettaras(Super(iK)).Idx(1,:)             = FirstI;
            
            for iZ = 1:length(khettaras(Super(iK)).P)
                % Also change all the individual values of the cells of the
                % khettaras.
                khettaras(Super(iK)).P(iZ).z   = [10^10 10^10];
                khettaras(Super(iK)).P(iZ).idx = FirstI;
            end
        end
    elseif S.SuperKhetFJ && S.SuperKhetHa && S.SuperKhetFezna
        Super = strmatchi('Super Khettara Fezna Jorf',{khettaras.name});
        Super = [Super strmatchi('Super Khettara Hannabou',{khettaras.name})];
        Super = [Super strmatchi('Super Khet Fezna',{khettaras.name})];
        for iK=1:length(khettaras)
            if any(iK == Super)
                continue;
            end
            HKe = strmatchi('H',khettaras(iK).vertexHdr);
            Z = strmatchi('Z',khettaras(iK).vertexHdr);

            % high values, so no water can flow into the drain (this is the
            % elevation),
            khettaras(iK).vertex(:,HKe) = 10^10;    
            khettaras(iK).vertex(:,Z(2)) = 10^10;
            khettaras(iK).cellLineVals(:,Z(2)) = 10^10;
            khettaras(iK).cellLineVals(:,HKe) = 10^10;
            khettaras(iK).Idx(1,:)             = FirstI;
            
            for iZ = 1:length(khettaras(iK).P)
                % Also change all the individual values of the cells of the
                % khettaras.
                khettaras(iK).P(iZ).z = [10^10 10^10];
                khettaras(iK).P(iZ).idx = FirstI;
            end
        end
    elseif S.SuperKhetFJ && S.SuperKhetHa
        Super = strmatchi('Super Khettara Fezna Jorf',{khettaras.name});
        Super = [Super strmatchi('Super Khettara Hannabou',{khettaras.name})];
        for iK=1:length(khettaras)
            if any(iK == Super)
                continue;
            end
            HKe = strmatchi('H',khettaras(iK).vertexHdr);
            Z = strmatchi('Z',khettaras(iK).vertexHdr);

            % high values, so no water can flow into the drain (this is the
            % elevation),
            khettaras(iK).vertex(:,HKe) = 10^10;    
            khettaras(iK).vertex(:,Z(2)) = 10^10;
            khettaras(iK).cellLineVals(:,Z(2)) = 10^10;
            khettaras(iK).cellLineVals(:,HKe) = 10^10;
            khettaras(iK).Idx(1,:)             = FirstI;
            
            for iZ = 1:length(khettaras(iK).P)
                % Also change all the individual values of the cells of the
                % khettaras.
                khettaras(iK).P(iZ).z = [10^10 10^10];
                khettaras(iK).P(iZ).idx = FirstI;
            end
        end
    elseif S.SuperKhetHa && S.SuperKhetFezna
        Super = strmatchi('Super Khettara Hannabou',{khettaras.name});
        Super = [Super strmatchi('Super Khet Fezna',{khettaras.name})];
        for iK=1:length(khettaras)
            if any(iK == Super)
                continue;
            end
            HKe = strmatchi('H',khettaras(iK).vertexHdr);
            Z = strmatchi('Z',khettaras(iK).vertexHdr);

            % high values, so no water can flow into the drain (this is the
            % elevation),
            khettaras(iK).vertex(:,HKe) = 10^10;    
            khettaras(iK).vertex(:,Z(2)) = 10^10;
            khettaras(iK).cellLineVals(:,Z(2)) = 10^10;
            khettaras(iK).cellLineVals(:,HKe) = 10^10;
            khettaras(iK).Idx(1,:)             = FirstI;
            
            for iZ = 1:length(khettaras(iK).P)
                % Also change all the individual values of the cells of the
                % khettaras.
                khettaras(iK).P(iZ).z = [10^10 10^10];
                khettaras(iK).P(iZ).idx = FirstI;
            end
        end
    elseif S.SuperKhetFJ && S.SuperKhetFezna
        Super = strmatchi('Super Khettara Fezna Jorf',{khettaras.name});
        Super = [Super strmatchi('Super Khet Fezna',{khettaras.name})];
        for iK=1:length(khettaras)
            if any(iK == Super)
                continue;
            end
            HKe = strmatchi('H',khettaras(iK).vertexHdr);
            Z = strmatchi('Z',khettaras(iK).vertexHdr);

            % high values, so no water can flow into the drain (this is the
            % elevation),
            khettaras(iK).vertex(:,HKe) = 10^10;    
            khettaras(iK).vertex(:,Z(2)) = 10^10;
            khettaras(iK).cellLineVals(:,Z(2)) = 10^10;
            khettaras(iK).cellLineVals(:,HKe) = 10^10;
            khettaras(iK).Idx(1,:)             = FirstI;
            
            for iZ = 1:length(khettaras(iK).P)
                % Also change all the individual values of the cells of the
                % khettaras.
                khettaras(iK).P(iZ).z = [10^10 10^10];
                khettaras(iK).P(iZ).idx = FirstI;
            end
        end
    elseif S.SuperKhetFJ
        Super = strmatchi('Super Khettara Fezna Jorf',{khettaras.name});
        for iK=1:length(khettaras)
            if any(iK == Super)
                continue;
            end
            HKe = strmatchi('H',khettaras(iK).vertexHdr);
            Z = strmatchi('Z',khettaras(iK).vertexHdr);

            % high values, so no water can flow into the drain (this is the
            % elevation),
            khettaras(iK).vertex(:,HKe) = 10^10;    
            khettaras(iK).vertex(:,Z(2)) = 10^10;
            khettaras(iK).cellLineVals(:,Z(2)) = 10^10;
            khettaras(iK).cellLineVals(:,HKe) = 10^10;
            khettaras(iK).Idx(1,:)             = FirstI;
            
            for iZ = 1:length(khettaras(iK).P)
                % Also change all the individual values of the cells of the
                % khettaras.
                khettaras(iK).P(iZ).z = [10^10 10^10];
                khettaras(iK).P(iZ).idx = FirstI;
            end
        end
    elseif S.SuperKhetHa
        Super = strmatchi('Super Khettara Hannabou',{khettaras.name});
        for iK=1:length(khettaras)
            if any(iK == Super)
                continue;
            end
            HKe = strmatchi('H',khettaras(iK).vertexHdr);
            Z = strmatchi('Z',khettaras(iK).vertexHdr);

            % high values, so no water can flow into the drain (this is the
            % elevation),
            khettaras(iK).vertex(:,HKe) = 10^10;    
            khettaras(iK).vertex(:,Z(2)) = 10^10;
            khettaras(iK).cellLineVals(:,Z(2)) = 10^10;
            khettaras(iK).cellLineVals(:,HKe) = 10^10;
            khettaras(iK).Idx(1,:)             = FirstI;
            
            for iZ = 1:length(khettaras(iK).P)
                % Also change all the individual values of the cells of the
                % khettaras.
                khettaras(iK).P(iZ).z = [10^10 10^10];
                khettaras(iK).P(iZ).idx = FirstI;
            end
        end
   elseif S.SuperKhetFezna % Theo 151014
        Super = strmatchi('Super',{khettaras.name});
        FeznaSuper = strmatchi('Super Khet Fezna',{khettaras.name});
        SuperOff = Super(~ismember(Super,FeznaSuper));
        for iK=SuperOff
            HKe = strmatchi('H',khettaras(iK).vertexHdr);
            Z = strmatchi('Z',khettaras(iK).vertexHdr);

            % high values, so no water can flow into the drain (this is the
            % elevation),
            khettaras(iK).vertex(:,HKe) = 10^10;    
            khettaras(iK).vertex(:,Z(2)) = 10^10;
            khettaras(iK).cellLineVals(:,Z(2)) = 10^10;
            khettaras(iK).cellLineVals(:,HKe) = 10^10;
            khettaras(iK).Idx(1,:)             = FirstI;
            
            for iZ = 1:length(khettaras(iK).P)
                % Also change all the individual values of the cells of the
                % khettaras.
                khettaras(iK).P(iZ).z = [10^10 10^10];
                khettaras(iK).P(iZ).idx = FirstI;
            end
        end
    elseif S.SuperKhetFezna && ~S.HoldKhettaras
        Super = strmatchi('Super Khet Fezna',{khettaras.name});
        for iK=1:length(khettaras)
            if any(iK == Super)
                continue;
            end
            HKe = strmatchi('H',khettaras(iK).vertexHdr);
            Z = strmatchi('Z',khettaras(iK).vertexHdr);

            % high values, so no water can flow into the drain (this is the
            % elevation),
            khettaras(iK).vertex(:,HKe) = 10^10;    
            khettaras(iK).vertex(:,Z(2)) = 10^10;
            khettaras(iK).cellLineVals(:,Z(2)) = 10^10;
            khettaras(iK).cellLineVals(:,HKe) = 10^10;
            khettaras(iK).Idx(1,:)             = FirstI;
            
            for iZ = 1:length(khettaras(iK).P)
                % Also change all the individual values of the cells of the
                % khettaras.
                khettaras(iK).P(iZ).z = [10^10 10^10];
                khettaras(iK).P(iZ).idx = FirstI;
            end
        end
    elseif S.SuperKhetFezna && S.HoldKhettaras
        Super = strmatchi('Super Khettara Fezna',{khettaras.name});
%        Super = [Super strmatchi('Super Khettara Hannabou',{khettaras.name})];
        for iK=1:length(khettaras)
            if any(iK == Super)
                HKe = strmatchi('H',khettaras(iK).vertexHdr);
                Z = strmatchi('Z',khettaras(iK).vertexHdr);

                % high values, so no water can flow into the drain (this is the
                % elevation),
                khettaras(iK).vertex(:,HKe) = 10^10;    
                khettaras(iK).vertex(:,Z(2)) = 10^10;
                khettaras(iK).cellLineVals(:,Z(2)) = 10^10;
                khettaras(iK).cellLineVals(:,HKe) = 10^10;
                khettaras(iK).Idx(1,:)             = FirstI;
                for iZ = 1:length(khettaras(iK).P)
                    % Also change all the individual values of the cells of the
                    % khettaras.
                    khettaras(iK).P(iZ).z = [10^10 10^10];
                    khettaras(iK).P(iZ).idx = FirstI;
                end               
            end
        end
    end
    
        FJKI = strmatchi('FeznaJorfKhet-Inf',{pumpAreas.name});
        HaKI = strmatchi('HannabouKhet-Inf',{pumpAreas.name});
        FKI  = strmatchi('FeznaKhet-Inf',{pumpAreas.name});

    % Here you change the layer of infiltration to the choosen layer in 
    % Koen_scenarios.m
    pumpAreas(FJKI) = changeLayer(pumpAreas(FJKI),S.LayerInf,gr,IBOUND);
    pumpAreas(HaKI) = changeLayer(pumpAreas(HaKI),S.LayerInf,gr,IBOUND);
        
        Flux = strmatchi('flux',pumpAreas.vertexHdr);
        Iflux= strmatchi('1flux',pumpAreas.vertexHdr);
        
        if exist('W','var') 
            % Using the outflow from the previous run (only if W 
            % (waterbalance) exist.) the values are - because it flows into
            % the model
            KRateFJ = - mean(W.QKhettaraFJ) * 1000000 / 365.24 / ...
                sum(pumpAreas(FJKI).A) * S.InfFact; 
            KRateHa = - mean(W.QKhettaraHa) * 1000000 / 365.24 / ...
                sum(pumpAreas(HaKI).A) * S.InfFact; 
            try
                KRateF  = - mean(W.QKhettaraF) * 1000000 / 365.24 / ...
                    sum(pumpAreas(FKI).A) * S.InfFact;
            catch
                KRateF = 0;
            end
        else
            KRateFJ = 0;
            KRateHa = 0;
            KRateF = 0;
        end
        
        % Here you set the infiltration rate to the areas.
        pumpAreas(FJKI).vertex(:,Flux) = KRateFJ;      
        pumpAreas(FJKI).vertex(:,Iflux) = KRateFJ;     
        pumpAreas(FJKI).cellLineVals(:,Flux) = KRateFJ;
        pumpAreas(FJKI).cellLineVals(:,Iflux) = KRateFJ;
        pumpAreas(FJKI).cellAreaVals(:,Flux) = KRateFJ;
        pumpAreas(FJKI).cellAreaVals(:,Iflux) = KRateFJ;
        
        pumpAreas(HaKI).vertex(:,Flux) = KRateHa;     
        pumpAreas(HaKI).vertex(:,Iflux) = KRateHa;     
        pumpAreas(HaKI).cellLineVals(:,Flux) = KRateHa;
        pumpAreas(HaKI).cellLineVals(:,Iflux) = KRateHa;
        pumpAreas(HaKI).cellAreaVals(:,Flux) = KRateHa;
        pumpAreas(HaKI).cellAreaVals(:,Iflux) = KRateHa;
        
        pumpAreas(FKI).vertex(:,Flux) = KRateF;     
        pumpAreas(FKI).vertex(:,Iflux) = KRateF;     
        pumpAreas(FKI).cellLineVals(:,Flux) = KRateF;
        pumpAreas(FKI).cellLineVals(:,Iflux) = KRateF;
        pumpAreas(FKI).cellAreaVals(:,Flux) = KRateF;
        pumpAreas(FKI).cellAreaVals(:,Iflux) = KRateF;
    
    % Only turned on if KhetLength = true. It will exclude parts that fall
    % outside of the S.Radius from the starting point of the khettara
    if S.KhetLength
        for iK = 1:length(khettaras)
            HKe = strmatchi('H',khettaras(iK).vertexHdr); % find H location
            Z   = strmatchi('Z',khettaras(iK).vertexHdr); % find Z location
            KX  = strmatchi('x',khettaras(iK).vertexHdr); % find X location
            KY  = strmatchi('y',khettaras(iK).vertexHdr); % find y location
            
            % for khettaras.vertex : 
            XKB = khettaras(iK).vertex(end,KX); % find starting point (x)
            YKB = khettaras(iK).vertex(end,KY); % find starting point (y)
            
            % check for every value the radius. if the radius is larger
            % then S.Radius, the values will be so high that no water can
            % infiltrate the drain.
            sizeK = size(khettaras(iK).vertex); % need y direction, because 
                                                % x can be larger than y.
            for iV = 1:sizeK(1)
                XK = khettaras(iK).vertex(iV,KX);
                YK = khettaras(iK).vertex(iV,KY);
                
                r = sqrt((XK - XKB)^2 - (YK - YKB)^2);
            
                if r > S.Radius
                    khettaras(iK).vertex(iV,HKe) = 10^10;
                    khettaras(iK).vertex(iV,Z(2)) = 10^10;
                end
            end
            
            % for khettaras.P & khettaras.cellLineVals :
            XKB = khettaras(iK).P(end).x(1,end); % Begin X of the Khettara
            YKB = khettaras(iK).P(end).y(1,end); % Begin Y of the Khettara

            for iP = 1:length(khettaras(iK).P)
                XK = khettaras(iK).P(iP).x(1,end); % X of khettara point
                YK = khettaras(iK).P(iP).y(1,end); % Y of khettara point

                r = sqrt((XK - XKB)^2 + (YK - YKB)^2); % radius from starting point

                % If the radius is longer than a certain value, the z value will be
                % so high the remaining khettara length will not drain any water
                if r > S.Radius
                    khettaras(iK).P(iP).z = [10^10 10^10];
                    khettaras(iK).cellLineVals(iP,Z(2)) = 10^10;
                    khettaras(iK).cellLineVals(iP,HKe) = 10^10;
                    
                    khettaras(iK).P(iP).idx = 1;
                    khettaras(iK).Idx(iP) = 1;
                end
            end
        end
    end
else
    % it will all be 0.
    khettaras = lineObj( gr,[XLSData, 'JorfData.xls'],'khettaras','DRN');
    pumpAreas = area2Obj(gr,[XLSData, 'JorfData.xls'],'Pumping','WEL');
    FirstI    = find(IBOUND,1); % Used to change Idx of khettaras to get 
                                % better plots
    
    if ~S.SuperKhetFJ && ~S.SuperKhetHa && ~S.SuperKhetFezna
        Super = strmatchi('Super',{khettaras.name});
        for iK=1:length(Super)
            HKe = strmatchi('H',khettaras(Super(iK)).vertexHdr);
            Z = strmatchi('Z',khettaras(Super(iK)).vertexHdr);
            C = strmatchi('c',khettaras(Super(iK)).vertexHdr);
            

            % high values, so no water can flow into the drain (this is the
            % elevation),
            khettaras(Super(iK)).vertex(:,HKe)        = 10^10;    
            khettaras(Super(iK)).vertex(:,Z(2))       = 10^10;
            khettaras(Super(iK)).cellLineVals(:,Z(2)) = 10^10;
            khettaras(Super(iK)).cellLineVals(:,HKe)  = 10^10;
            khettaras(Super(iK)).vertex(:,Z(1))       = 10^10;
            khettaras(Super(iK)).cellLineVals(:,Z(1)) = 10^10;
            khettaras(Super(iK)).vertex(:,C(1))       = 10^10; 
            khettaras(Super(iK)).cellLineVals(:,C(1)) = 10^10;
            khettaras(Super(iK)).Idx(1,:)             = FirstI;
            
            for iZ = 1:length(khettaras(Super(iK)).P)
                % Also change all the individual values of the cells of the
                % khettaras.
                khettaras(Super(iK)).P(iZ).z   = [10^10 10^10];
                khettaras(Super(iK)).P(iZ).idx = FirstI;
            end
        end
    elseif S.SuperKhetFJ && S.SuperKhetHa && S.SuperKhetFezna
        Super = strmatchi('Super Khettara Fezna Jorf',{khettaras.name});
        Super = [Super strmatchi('Super Khettara Hannabou',{khettaras.name})];
        Super = [Super strmatchi('Super Khet Fezna',{khettaras.name})];
        for iK=1:length(khettaras)
            if any(iK == Super)
                continue;
            end
            HKe = strmatchi('H',khettaras(iK).vertexHdr);
            Z = strmatchi('Z',khettaras(iK).vertexHdr);

            % high values, so no water can flow into the drain (this is the
            % elevation),
            khettaras(iK).vertex(:,HKe) = 10^10;    
            khettaras(iK).vertex(:,Z(2)) = 10^10;
            khettaras(iK).cellLineVals(:,Z(2)) = 10^10;
            khettaras(iK).cellLineVals(:,HKe) = 10^10;
            khettaras(iK).Idx(1,:)             = FirstI;
            
            for iZ = 1:length(khettaras(iK).P)
                % Also change all the individual values of the cells of the
                % khettaras.
                khettaras(iK).P(iZ).z = [10^10 10^10];
                khettaras(iK).P(iZ).idx = FirstI;
            end
        end
    elseif S.SuperKhetFJ && S.SuperKhetHa
        Super = strmatchi('Super Khettara Fezna Jorf',{khettaras.name});
        Super = [Super strmatchi('Super Khettara Hannabou',{khettaras.name})];
        for iK=1:length(khettaras)
            if any(iK == Super)
                continue;
            end
            HKe = strmatchi('H',khettaras(iK).vertexHdr);
            Z = strmatchi('Z',khettaras(iK).vertexHdr);

            % high values, so no water can flow into the drain (this is the
            % elevation),
            khettaras(iK).vertex(:,HKe) = 10^10;    
            khettaras(iK).vertex(:,Z(2)) = 10^10;
            khettaras(iK).cellLineVals(:,Z(2)) = 10^10;
            khettaras(iK).cellLineVals(:,HKe) = 10^10;
            khettaras(iK).Idx(1,:)             = FirstI;
            
            for iZ = 1:length(khettaras(iK).P)
                % Also change all the individual values of the cells of the
                % khettaras.
                khettaras(iK).P(iZ).z = [10^10 10^10];
                khettaras(iK).P(iZ).idx = FirstI;
            end
        end
    elseif S.SuperKhetHa && S.SuperKhetFezna
        Super = strmatchi('Super Khettara Hannabou',{khettaras.name});
        Super = [Super strmatchi('Super Khet Fezna',{khettaras.name})];
        for iK=1:length(khettaras)
            if any(iK == Super)
                continue;
            end
            HKe = strmatchi('H',khettaras(iK).vertexHdr);
            Z = strmatchi('Z',khettaras(iK).vertexHdr);

            % high values, so no water can flow into the drain (this is the
            % elevation),
            khettaras(iK).vertex(:,HKe) = 10^10;    
            khettaras(iK).vertex(:,Z(2)) = 10^10;
            khettaras(iK).cellLineVals(:,Z(2)) = 10^10;
            khettaras(iK).cellLineVals(:,HKe) = 10^10;
            khettaras(iK).Idx(1,:)             = FirstI;
            
            for iZ = 1:length(khettaras(iK).P)
                % Also change all the individual values of the cells of the
                % khettaras.
                khettaras(iK).P(iZ).z = [10^10 10^10];
                khettaras(iK).P(iZ).idx = FirstI;
            end
        end
    elseif S.SuperKhetFJ && S.SuperKhetFezna
        Super = strmatchi('Super Khettara Fezna Jorf',{khettaras.name});
        Super = [Super strmatchi('Super Khet Fezna',{khettaras.name})];
        for iK=1:length(khettaras)
            if any(iK == Super)
                continue;
            end
            HKe = strmatchi('H',khettaras(iK).vertexHdr);
            Z = strmatchi('Z',khettaras(iK).vertexHdr);

            % high values, so no water can flow into the drain (this is the
            % elevation),
            khettaras(iK).vertex(:,HKe) = 10^10;    
            khettaras(iK).vertex(:,Z(2)) = 10^10;
            khettaras(iK).cellLineVals(:,Z(2)) = 10^10;
            khettaras(iK).cellLineVals(:,HKe) = 10^10;
            khettaras(iK).Idx(1,:)             = FirstI;
            
            for iZ = 1:length(khettaras(iK).P)
                % Also change all the individual values of the cells of the
                % khettaras.
                khettaras(iK).P(iZ).z = [10^10 10^10];
                khettaras(iK).P(iZ).idx = FirstI;
            end
        end
    elseif S.SuperKhetFJ
        Super = strmatchi('Super Khettara Fezna Jorf',{khettaras.name});
        for iK=1:length(khettaras)
            if any(iK == Super)
                continue;
            end
            HKe = strmatchi('H',khettaras(iK).vertexHdr);
            Z = strmatchi('Z',khettaras(iK).vertexHdr);

            % high values, so no water can flow into the drain (this is the
            % elevation),
            khettaras(iK).vertex(:,HKe) = 10^10;    
            khettaras(iK).vertex(:,Z(2)) = 10^10;
            khettaras(iK).cellLineVals(:,Z(2)) = 10^10;
            khettaras(iK).cellLineVals(:,HKe) = 10^10;
            khettaras(iK).Idx(1,:)             = FirstI;
            
            for iZ = 1:length(khettaras(iK).P)
                % Also change all the individual values of the cells of the
                % khettaras.
                khettaras(iK).P(iZ).z = [10^10 10^10];
                khettaras(iK).P(iZ).idx = FirstI;
            end
        end
    elseif S.SuperKhetHa
        Super = strmatchi('Super Khettara Hannabou',{khettaras.name});
        for iK=1:length(khettaras)
            if any(iK == Super)
                continue;
            end
            HKe = strmatchi('H',khettaras(iK).vertexHdr);
            Z = strmatchi('Z',khettaras(iK).vertexHdr);

            % high values, so no water can flow into the drain (this is the
            % elevation),
            khettaras(iK).vertex(:,HKe) = 10^10;    
            khettaras(iK).vertex(:,Z(2)) = 10^10;
            khettaras(iK).cellLineVals(:,Z(2)) = 10^10;
            khettaras(iK).cellLineVals(:,HKe) = 10^10;
            khettaras(iK).Idx(1,:)             = FirstI;
            
            for iZ = 1:length(khettaras(iK).P)
                % Also change all the individual values of the cells of the
                % khettaras.
                khettaras(iK).P(iZ).z = [10^10 10^10];
                khettaras(iK).P(iZ).idx = FirstI;
            end
        end
    elseif S.SuperKhetFezna && ~S.HoldKhettaras
        Super = strmatchi('Super Khet Fezna',{khettaras.name});
        for iK=1:length(khettaras)
            if any(iK == Super)
                continue;
            end
            HKe = strmatchi('H',khettaras(iK).vertexHdr);
            Z = strmatchi('Z',khettaras(iK).vertexHdr);

            % high values, so no water can flow into the drain (this is the
            % elevation),
            khettaras(iK).vertex(:,HKe) = 10^10;    
            khettaras(iK).vertex(:,Z(2)) = 10^10;
            khettaras(iK).cellLineVals(:,Z(2)) = 10^10;
            khettaras(iK).cellLineVals(:,HKe) = 10^10;
            khettaras(iK).Idx(1,:)             = FirstI;
            
            for iZ = 1:length(khettaras(iK).P)
                % Also change all the individual values of the cells of the
                % khettaras.
                khettaras(iK).P(iZ).z = [10^10 10^10];
                khettaras(iK).P(iZ).idx = FirstI;
            end
        end
    elseif S.SuperKhetFezna && S.HoldKhettaras
        Super = strmatchi('Super Khettara Fezna Jorf',{khettaras.name});
        Super = [Super strmatchi('Super Khettara Hannabou',{khettaras.name})];
        for iK=1:length(khettaras)
            if any(iK == Super)
                HKe = strmatchi('H',khettaras(iK).vertexHdr);
                Z = strmatchi('Z',khettaras(iK).vertexHdr);

                % high values, so no water can flow into the drain (this is the
                % elevation),
                khettaras(iK).vertex(:,HKe) = 10^10;    
                khettaras(iK).vertex(:,Z(2)) = 10^10;
                khettaras(iK).cellLineVals(:,Z(2)) = 10^10;
                khettaras(iK).cellLineVals(:,HKe) = 10^10;
                khettaras(iK).Idx(1,:)             = FirstI;
                for iZ = 1:length(khettaras(iK).P)
                    % Also change all the individual values of the cells of the
                    % khettaras.
                    khettaras(iK).P(iZ).z = [10^10 10^10];
                    khettaras(iK).P(iZ).idx = FirstI;
                end               
            end
        end
    end
    
    for iK = 1:length(khettaras)
        HKe = strmatchi('H',khettaras(iK).vertexHdr);
        Z = strmatchi('Z',khettaras(iK).vertexHdr);
        C = strmatchi('c',khettaras(iK).vertexHdr);
        
        % high values, so no water can flow into the drain (this is the
        % elevation),
        khettaras(iK).vertex(:,HKe) = 10^10;    
        khettaras(iK).vertex(:,Z(2)) = 10^10;
        khettaras(iK).cellLineVals(:,Z(2)) = 10^10;
        khettaras(iK).cellLineVals(:,HKe) = 10^10;
        khettaras(iK).vertex(:,C(1)) = 10^10; 
        khettaras(iK).cellLineVals(:,C(1)) = 10^10;
        
        for iZ = 1:length(khettaras(iK).P)
            % Also change all the individual values of the cells of the
            % khettaras.
            khettaras(iK).P(iZ).z = [10^10 10^10];
        end
    end
    
    Flux = strmatchi('flux',pumpAreas.vertexHdr);
    Iflux= strmatchi('1flux',pumpAreas.vertexHdr);   
        
    FJKI = strmatchi('FeznaJorfKhet-Inf',{pumpAreas.name});
        pumpAreas(FJKI).vertex(:,Flux) = 0;     
        pumpAreas(FJKI).vertex(:,Iflux) = 0;     
        pumpAreas(FJKI).cellLineVals(:,Flux) = 0;
        pumpAreas(FJKI).cellLineVals(:,Iflux) = 0;
        pumpAreas(FJKI).cellAreaVals(:,Flux) = 0;
        pumpAreas(FJKI).cellAreaVals(:,Iflux) = 0;
    
    HaKI = strmatchi('HannabouKhet-Inf',{pumpAreas.name});
        pumpAreas(HaKI).vertex(:,Flux) = 0;      
        pumpAreas(HaKI).vertex(:,Iflux) = 0;     
        pumpAreas(HaKI).cellLineVals(:,Flux) = 0;
        pumpAreas(HaKI).cellLineVals(:,Iflux) = 0;
        pumpAreas(HaKI).cellAreaVals(:,Flux) = 0;
        pumpAreas(HaKI).cellAreaVals(:,Iflux) = 0;
        
    FKI  = strmatchi('FeznaKhet-Inf',{pumpAreas.name});
        pumpAreas(FKI).vertex(:,Flux) = 0;     
        pumpAreas(FKI).vertex(:,Iflux) = 0;     
        pumpAreas(FKI).cellLineVals(:,Flux) = 0;
        pumpAreas(FKI).cellLineVals(:,Iflux) = 0;
        pumpAreas(FKI).cellAreaVals(:,Flux) = 0;
        pumpAreas(FKI).cellAreaVals(:,Iflux) = 0;
end

if 1 
    % Always true, pumping is defined in scenarios. If no pumping Pumprate 
    % will be 0;
    Flux = strmatchi('flux',pumpAreas.vertexHdr);
    Iflux= strmatchi('1flux',pumpAreas.vertexHdr);  
        
    FJP  = strmatchi('FeznaJorfPump',{pumpAreas.name});
    FJPI = strmatchi('FeznaJorfPum-Inf',{pumpAreas.name});
    OHP  = strmatchi('Oukhit',{pumpAreas.name});
    OHPI = strmatchi('Oukhi-inf',{pumpAreas.name});
    
    if S.Pumping
        % defined in scenarios
        PumpRate = -(S.PumpVol) / 365.24 / sum(pumpAreas(FJP).A); 
    else
        PumpRate = 0;
    end
    
    % change the layer of infiltration.
    pumpAreas(FJPI) = changeLayer(pumpAreas(FJPI),S.LayerInf,gr,IBOUND);
         pumpAreas(FJP).vertex(:,Flux) = PumpRate;       
         pumpAreas(FJP).vertex(:,Iflux) = PumpRate;
         pumpAreas(FJP).cellLineVals(:,Flux) = PumpRate; 
         pumpAreas(FJP).cellLineVals(:,Iflux) = PumpRate;
         pumpAreas(FJP).cellAreaVals(:,Flux) = PumpRate; 
         pumpAreas(FJP).cellAreaVals(:,Iflux) = PumpRate;
         
    if S.Oukhit
        PumpRateOH = -(S.PumpVolOuk) / 365.24 / sum(pumpAreas(OHP).A); 
    else
        PumpRateOH = 0;
    end
    
    % change the layer of infiltration.
    pumpAreas(OHPI) = changeLayer(pumpAreas(OHPI),S.LayerInf,gr,IBOUND);
         pumpAreas(OHP).vertex(:,Flux) = PumpRateOH;       
         pumpAreas(OHP).vertex(:,Iflux) = PumpRateOH;
         pumpAreas(OHP).cellLineVals(:,Flux) = PumpRateOH; 
         pumpAreas(OHP).cellLineVals(:,Iflux) = PumpRateOH;
         pumpAreas(OHP).cellAreaVals(:,Flux) = PumpRateOH; 
         pumpAreas(OHP).cellAreaVals(:,Iflux) = PumpRateOH;

    if exist('W','var') 
        % Using the outflow from the previous run (only if W (waterbalance)
        % exist) it is - because it flows into the model
        PumpInf = - mean(W.QFJP) * 1000000 / 365.24 / ...
            sum(pumpAreas(FJPI).A) * S.InfFact;
        try 
            % Because Oukhit has been added later, it is possible it is not 
            % in the previous saved waterbalance in the storage file
            PumpInfOH = - mean(W.QOHP) * 1000000 / 365.24 / ...
                sum(pumpAreas(OHPI).A) * S.InfFact;
        catch
            PumpInfOH = 0;
        end
    elseif S.LayerInf == 3.5
        % if W does not exist and infiltration is in the deepest layer, the
        % rate needs to be defined so it can be used in the waterbalance,
        % in the next run the PumpInf can be calculated with the
        % information of the waterbalance. By doing this it will get to
        % equilibrium faster.
        PumpInf = PumpRate * S.InfFact;
    else
        PumpInf   = 0;
        PumpInfOH = 0;
    end
        pumpAreas(FJPI).vertex(:,Flux) = PumpInf;       
        pumpAreas(FJPI).vertex(:,Iflux) = PumpInf;
        pumpAreas(FJPI).cellLineVals(:,Flux) = PumpInf; 
        pumpAreas(FJPI).cellLineVals(:,Iflux) = PumpInf;   
        pumpAreas(FJPI).cellAreaVals(:,Flux) = PumpInf; 
        pumpAreas(FJPI).cellAreaVals(:,Iflux) = PumpInf; 
        
        pumpAreas(OHPI).vertex(:,Flux) = PumpInfOH;       
        pumpAreas(OHPI).vertex(:,Iflux) = PumpInfOH;
        pumpAreas(OHPI).cellLineVals(:,Flux) = PumpInfOH; 
        pumpAreas(OHPI).cellLineVals(:,Iflux) = PumpInfOH;   
        pumpAreas(OHPI).cellAreaVals(:,Flux) = PumpInfOH; 
        pumpAreas(OHPI).cellAreaVals(:,Iflux) = PumpInfOH;  
end

if 1 
    % allways true, Flooding is defined in scenarios. S.FloodingFJ will be 
    % 0 is no flooding is selected
    FloodAreas = area2Obj(gr,[XLSData, 'JorfData.xls'],'Floods','WEL');
    FJFI = strmatchi('FeznaJorfFlood-inf',{FloodAreas.name});
    HaFI = strmatchi('Hannabou-Flood',{FloodAreas.name});
    Sand = strmatchi('SandInf',{FloodAreas.name});
    
%     FloodAreas(FJFI) = changeLayer(FloodAreas(FJFI),S.LayerInf,gr,IBOUND);
%     FLoodAreas(HaFI) = changeLayer(FloodAreas(HaFI),S.LayerInf,gr,IBOUND);
    
    if S.Floods
        % here you multiply the recharge with a constant, in order to get right
        % amount of flooding. First the constant is calculated, then it is
        % used to get the correct factor.
        
        FJFloodFactor = (S.FloodingFJ / sum(FloodAreas(FJFI).A) / 365.24 ) /...
            sum(S.RainRech);
        HaFloodFactor = (S.FloodingHa / sum(FloodAreas(HaFI).A) / 365.24 ) /...
            sum(S.RainRech);
        SandFloodFactor = (S.FloodingSand / (sum(FloodAreas(Sand(1)).A) + ...
            sum(FloodAreas(Sand(2)).A) + sum(FloodAreas(Sand(3)).A)) ...
            /365.24 )/sum(S.RainRech);
        
        facFlood = (S.PercentFlooding * ones(1,12)).';
        facFlood = mean(facFlood(:)); % for the time being not yet split over all years
        FloodAreas(FJFI).V = num2cell(FJFloodFactor .* facFlood(:) ...
             * S.InfFact / perlen(1));
        FloodAreas(HaFI).V = num2cell(HaFloodFactor .* facFlood(:) ...
             * S.InfFact / perlen(1));
        
        for iS = 1:length(Sand)
            if S.SandInf
                FloodAreas(Sand(iS)).V = num2cell(SandFloodFactor ...
                    .* facFlood(:) * S.InfFact / perlen(1));
            else
                FloodAreas(Sand(iS)).V = {0};
            end
        end
    else
        for iF = 1:length(FloodAreas)
            FloodAreas(iF).V = {0};
        end
    end
end

if S.Piezom
    piezom    = pointObj(gr,[XLSData,'JorfData.xls'],'piezom','NIL');
end



%% LAYER THICKNESS
% Show plots of the model so far.
if S.LayerThick
    D = sum(gr.DZ,3);

    figure('name',sprintf(' %s',esriLayers{:}),'pos',screenPos(0.75));
    axes('nextplot','add');
    title(sprintf('total thickness of model <<%s>>',sprintfs(' + %s',esriLayers(2:end))))  ;
    xlabel('x [m]'); ylabel('y [m]');

    [c,h] = gr.contourf(D .* ibound,0:1:50);
    clabel(c,h,'color','w');
    h=colorbar; set(get(h,'title'),'string','Total thickness [m]');

    set(gcf, 'renderer', 'zbuffer');
end

%% TOTAL TRANSMISSIVITY

if S.TotalTrans
    kD = gr.DZ .* HK;

    figure('name','total transmissivity [m2/d]','pos',screenPos(0.75));
    axes('nextplot','add');
    title(sprintf('total kD of model <<%s>>',sprintfs(' + %s',esriLayers(2:end))));
    xlabel('x [m]'); ylabel('y [m]');
    kDrange = contRange(kD,25);

    [c,h] = gr.contourf(sum(kD,3),kDrange); clabel(c,h,'color','y');
    h=colorbar; set(get(h,'title'),'string','Total kD [m2/d]');
    set(gcf, 'renderer', 'zbuffer');
end

%% Save pass to mf_analyze

save underneath XLSData catchment

fprintf('***** %s -- finished. Running it took %g seconds *****\n',mfilename,toc);
