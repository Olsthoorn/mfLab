function well = concXY(o,varargin)
%ANIMATEOBJ/CONCXY -- plot concentrations in the horizontal plane
%
% USAGE:
%    well       = animateObj.concXZ([gr,],speciesIndex|speciesNames,ILay,well,...
%                   options,value,option,value,option,value);
%
% gr       = gridObj
% species  = index/Nr of compoments/Species to plot.
%            Names of components/species may be used instead, but must be given as
%            a cell array like {'tracer','temperature'}. These species names
%            must be in the order used in the model.
%            If the concentrations should be overlain by head contours,
%            add 'head' at the end like this: {'tracer',temperature','head'}
%            These species names must have been specied when calling animateObj
% Ilay     = ILay is the layer range whose concentration will be averaged
%            before plotting. ILay is all if omitted.
%            ILay may also be a logical expression.
% well       are of class wellsObj or of class MNE1Obj to be included in the plot
%
% varargin can be used to specify further options
%
% options that can be used as option,value pairs. Options can be used in
% stead of the call under USAGE. If options are used, their position in the
% call is immaterial.
%    'species' Numbers of the species to be plotted or 
%              {cell array with names of species}
%               add 'head' to plot head contours over the concentrations.
%    'contourClr' clr    to specify the color of the head contours
%    'ILay  ',layer numbers to be averaged vertically
%    'figure' figname
%    'figPos' figpos [xLL yLL w h] in pixels
%    'STCONC' STCONC  % start conc, STCONC{n} must contain the values for
%             species n
%    'xlim' xlim
%    'zlim' zlim (vertical dimension to show)
%    'VK',VK to compute and visualize the hydraulic resistance of layers.
%    'backgr', mfile  with plotting instructions to generate lines etc on top
%               of the picture being drawn
%
%  varargin can also be used to pass plotting property,value pairs
%  it depends on the receiving Matlab function which work.
%
% EXAMPLES
%   well = animateObj.concXZ(gr,{'Salinity','Temp'},ILay,well,varargin);
%
% TO 120908 121017 130403 130616

aYear = 365.24; % days per year

%% Make sure there is a grid in the input or there is a grid in animateObj
%  This allows specifying in this function and overwriting in animateObj.
%  In the end there must be a gridObj.
[o.gr,varargin] = getProp(varargin,'gridObj',[]);
if isempty(o.gr)
    [o.gr,varargin] = getType(varargin,'gridObj',[]);
end

if isempty(o.gr)
    error('%s: gridObj missing, is required because animateObj does not yet have one',mfilename);
end

%% First deal with all property value pairs in the input

% STCONC given? This is seen as request for conc difference from start values
[STCONC ,varargin] = getProp(varargin,'STCONC',[]);

Delta = ~isempty(STCONC);
if Delta && ~iscell(STCONC)
    prefix='Delta';
    STCONC = {STCONC};
else
    prefix='';
end

% New fig?
[figName,varargin] = getProp(varargin,'fig','');
[figPos ,varargin] = getProp(varargin,{'figPos','pos'},screenPos(0.75));

% plot the grid?
[gridLineSpec,varargin] = getProp(varargin,'plotgrid','c');

[xlim   ,varargin] = getProp(varargin,'xlim',o.gr.xm([1 end]));
[ylim   ,varargin] = getProp(varargin,'ylim',o.gr.ym([end,1]));

[xLbl   ,varargin] = getProp(varargin,'xlabel','x [m]');
[yLbl   ,varargin] = getProp(varargin,'ylabel','y [m]');
[zLbl   ,varargin] = getProp(varargin,'zlabel','z [m]');

% Background, through drawing instructions in script backgr.m
% These will be executed at te and of the plotting in time step 1, so that
% the results well be visual on top of the concentrations.
% Make sure backg.m works independently. Test it in the workspace first.
[backgr,varargin] = getProp(varargin,'backg',[]);

% A map from google maps requested (needs testing)
[map    ,varargin] = getProp(varargin,'map','roadmap'); % gooogle Earth picture type

% What coordinates, we have only 'rd' (Dutch) and 'wgs' (E,N or Lat,Lon)
[coords ,varargin] = getProp(varargin,'coord','rd');    % coordinate system 'rd or wgs

% color of head contours lines
[contourClr,varargin] = getProp(varargin,'contourClr',grey);

%% What is left-over should be clean input
[species,varargin] = getProp(varargin,'species',{});
if isempty(species)
    [species,varargin] = getNext(varargin,'cell',{});
end
if ischar(species), species={species}; end

if ~isempty(species)
    if iscell(species)
        IComp = strmatchi(species,o.titleStr);
        if IComp==0
            error('species names {%s} does not match species in animateObj {%s}',...
                sprintf(' %s',species{:}),sprintf(' %s',o.parent.titleStr{:})); %#ok
        end
    else % species was numeric
        IComp = species;
        IComp(IComp>numel(o.titleStr))=numel(o.titleStr);
        IComp(IComp<1)=1;
        IComp = unique(IComp);
    end
else
    [IComp  ,varargin] = getNext(varargin,{'double','logical'},1);
    if islogical(IComp)
        IComp = find(IComp);
    end
    IComp(IComp>numel(o.titleStr))=numel(o.titleStr);
    IComp(IComp<1) = 1;
    IComp = unique(IComp);
end

if isempty(IComp) || IComp==0;
    if o.NCOMP==1
        IComp=1;
    else
        error('%s: IComp not specified in call to %s',mfilename,mfilename);
    end
end

%% Specify a specific layer number?
[ILay     ,varargin] = getProp(varargin,'Layer',[]);
if ~isempty(ILay)
    if any(ILay<o.gr.zGr(end) || ILay>o.gr.zGr(1))
        error('%s: ILay must be between %d and %d.',mfilename,1,gr.Nlay);
    end
else   
    [ILay   ,varargin] = getNext(varargin,'double',[]);
end

if isempty(ILay)
    [ILay   ,varargin] = getNext(varargin,'logical',ILay);
    if ~isempty(ILay), ILay= find(ILay); end
end

if isempty(ILay)
    ILay=1;
end

if numel(size(ILay))>2
    error(['%s: ILay must be a vector, perhaps to forgot to specify it in call to %s\n',...
        'Probable cause: read in STCONC because ILay or IComp were not specified'],...
        mfilename,mfilename);
end

%% wells
[well,varargin]     = getType(varargin,{'wellObj','MNW1Obj','MNW2Obj'},[]);


%% assert presence of UCN
if isempty(o.UCN)
    error(['%s: No concenration data loaded to contour.\n',...
        'REMEDY: request animateObj with title of species corresponding to\n',...
        '        MT3DMS/SEAWAT output concentratioin files UNC000?.UCN.'],...
        mfilename);
else
    if ~all(o.gr.size == [size(o.UCN{1}(1).values,1),size(o.UCN{1}(1).values,2),size(o.UCN{1}(1).values,3)])
        error('size of grid [%d %d %d] does not match that of concentration [%d %d %d]',...
            mfilename,o.gr.size,size(o.UCN{1}.values));
    end
end

%% Verify size of STCONC
if Delta   
    if ~all(o.gr.size  == [size(STCONC{1},1),size(STCONC{1},2),size(STCONC{1},3)])
       error('%s: size of arrays in STCONC [%d %d %d] does not match that grid [%d %d %d]  in animate',...
           mfilename,size(STCONC{1}),o.gr.size);
    end    
end

%% assert presence of requrested UCN
if max(IComp) > numel(o.UCN)
    error('%s: requested compoment number <<%d>> is > number of components in aniamteObj <<%d>>',...
        mfilename,max(IComp),numel(o.UCN));
end

%% New figure explicitly requested
if ~isempty(figName)
    figure('name',figName,'position',figPos);
else
    if isempty(get(0,'children'))
        figure('name',figName,'position',figPos);
    else
        set(gcf,'position',figPos);
    end
end

%% default axis properties
axProps = {'nextplot','add','xgrid','on','ygrid','on','zgrid','on','xlim',xlim,'ylim',ylim};

crange{max(IComp)}=NaN;
clim{max(IComp)}  =NaN;
ax(max(IComp))    =NaN;

for iComp = IComp(:)'
    
    if Delta
        crange{iComp} = ContourRange({o.UCN{iComp},STCONC{iComp}},o.numberOfConcContours);
        clim{iComp} = crange{iComp}([1 end]);
    else
        crange{iComp} = ContourRange(o.UCN{iComp},o.numberOfConcContours);
        clim{iComp} = [0 crange{iComp}(end)];
    end

    ax(iComp) = subplot(o.NCOMP,1,iComp,axProps{:},'color','none','clim',clim{iComp},varargin{:});

    xlabel(ax(iComp),xLbl);
    ylabel(ax(iComp),yLbl);
    zlabel(ax(iComp),zLbl);
       
    %% Overlay google image
    try
        googleMap(ax(iComp),xlim,ylim,o.gr.zGr(end),map,coords,'rgb');
    catch %#ok
        warning('mfLab:concXY:googleMapLoadFailure',...
            '%s: Couldn''t get google Map, ... continue wihtout ...',mfilename);
    end
    
    axis('equal');
    axis('tight');
    %view(3);

    %% Plot the grid ??
    if ~isempty(gridLineSpec)
        if isLineSpec(gridLineSpec)
            o.gr.plotGrid(ax(iComp),gridLineSpec);
        end
    end

end

%% Set up movie
myName = mfilename; myName(1) = upper(myName(1));

vidObj = VideoWriter([o.basename myName]);
vidObj.FrameRate = o.framerate;
vidObj.Quality   = o.quality;
vidObj.open();

%% make movie
ht = NaN(o.NCOMP,1);
hC = NaN(o.NCOMP,1);

if ~isempty(well), wells{o.NCOMP}(numel(well),1)=wellObj(); end

dateformat = 'yyyy-mm-dd';
date0 = datestr(o.t0,dateformat);

time = o.t0+o.time;  % also works when time starts at 0

for it=1:length(time);
    
    for iComp=1:o.NCOMP
        if o.t0>0
            % if o.t>0 assume o.t is starting date as datenum
            tts = [prefix  o.titleStr{iComp} ',     ', date0, ' - ', datestr(time(it),dateformat)];
        else
            % if o.t==0 then use time as days
            if time(it)< 5 * aYear
                tts = [prefix  o.titleStr{iComp} ',     ', sprintf('t = %.4g d',time(it))];
            else
                tts = [prefix  o.titleStr{iComp} ',     ', sprintf('t = %.4g y',time(it)/aYear)];
            end
        end
        
        if it==1
            
            ht(iComp) = title(ax(iComp),tts);

            % concentration change
            if Delta
                [~,hC(iComp)] = contourf(ax(iComp),o.gr.xc,o.gr.yc,...
                    mean(o.UCN{iComp}(it).values(:,:,ILay)-STCONC{iComp}(:,:,ILay),3),...
                       crange{iComp},'edgecolor','none');
                   
                   set(ax(iComp),'clim',crange{iComp}([1 end]));
           % concentration
            else
                [~,hC(iComp)] = contourf(ax(iComp),o.gr.xc,o.gr.yc,...
                mean(o.UCN{iComp}(it).values(:,:,ILay),3),...
                   crange{iComp},'edgecolor','none');
            end
            
            % color bar
            hb = colorbar('peer',ax(iComp)); set(get(hb,'title'),'string',o.titleStr(iComp));
                  
            % if head is requested, draw contour lines
            if ~isempty(species) && strmatchi('head',species)
                [~,hh] = contour(ax(IComp),o.gr.xm,o.gr.ym,mean(o.H(it).values(:,:,ILay),3),o.hrange,contourClr);
            end

            % different wells because of different axes with their own handles
            if ~isempty(well)
                wells{iComp} = well.plot3D(ax(iComp),'color','w');
                wells{iComp} = wells{iComp}.plotXY(ax(iComp),'color','w');
            end            
            
            % background image?
            if ~isempty(backgr)
                    eval(backgr);
            end

        else
            set(ht(iComp),'string',tts);
            
            % update concentration change
            if Delta
                set(hC(iComp),'zdata' ,mean(o.UCN{iComp}(it).values(:,:,ILay)-STCONC{iComp}(:,:,ILay),3));
            % update concentration
            else
                set(hC(iComp),'zdata' ,mean(o.UCN{iComp}(it).values(:,:,ILay),3));
            end
            
            % update head contours
            if strmatchi('head',species)
                set(hh,'zdata',mean(o.H(it).values(:,:,ILay),3));
            end
            
            % update wells
           if ~isempty(well)
               wells{iComp} = wells{iComp}.plotXY(it);
           end
        end
    end
            
    vidObj.writeVideo(getframe(gcf));

end

vidObj.close();
