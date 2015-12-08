function o = slices(o,varargin)
 % aninmateObj.bubles(o.gr,IComp,well,STCONC,sx,sy,sz)
 % animate bubbles if STCONC is given, subtract them first
 %
 % TO 120914 130406

[figName,varargin] = getProp(varargin,'fig','');
[figPos ,varargin] = getProp(varargin,{'figpos','pos'},screenPos(0.75));
[ax     ,varargin] = getProp(varargin,'ax', []);
[map    ,varargin] = getProp(varargin,'map','roadmap');
[coords ,varargin] = getProp(varargin,'co','rd');
[grdLineSpec,varargin] = getProp(varargin,'plotGrid',[]);

if ~isempty(o.gr); % ignore gr in input
    [~,varargin] = getNext(varargin,'gridObj',[]);
else   % then get the grid
    [o.gr,varargin] = getNext(varargin,'gridObj',[]);
    if isempty(o.gr)
        error('%s: gridObj missing, is required because animateObj does not yet have one',mfilename);
    end
end

[IComp   ,varargin] = getNext(varargin,'double',1:o.NCOMP);

[well    ,varargin] = getNext(varargin,'wellObj',[]);
[well    ,varargin] = getNext(varargin,'MNW1Obj',well);
[well    ,varargin] = getNext(varargin,'MNW2Obj',well);

[STCONC ,varargin] = getNext(varargin,'numeric',[]);
[STCONC ,varargin] = getNext(varargin,'cell',STCONC);
[STCONC ,varargin] = getProp(varargin,'STCONC',STCONC);

if ~isempty(STCONC);
    Delta = true;
    prefix = 'Delta';
    if ~iscell(STCONC)
        STCONC = {STCONC};
    end
else
    Delta=false;
    prefix='';
end

[xlim, varargin] = getProp(varargin,'xlim',o.gr.xGr([1 end]));
[ylim, varargin] = getProp(varargin,'ylim',o.gr.yGr([end 1]));
[zlim, varargin] = getProp(varargin,'zlim',o.gr.zGr([end 1]));

[sx, varargin] = getNext(varargin,'double',[]);
[sy, varargin] = getNext(varargin,'double',[]);
[sz, varargin] = getNext(varargin,'double',[]);

if any(sx<o.gr.xGr(1) | sx>o.gr.xGr(end))
    error('%s: sx [%s ] not all within coordinates of x-axis [%g %g]',...
        mfilename,sprintf(' %g',sx),o.gr.xGr([1 end]));
end
if any(sy<o.gr.yGr(end) | sy>o.gr.yGr(1))
    error('%s: sy [%s ] not all within coordinates of y-axis [%g %g]',...
        mfilename,sprintf(' %g',sy),o.gr.yGr([end 1]));
end
if any(sz<o.gr.zGr(end) | sz>o.gr.zGr(1))
    error('%s: sz [%s ] not all within coordinates of z-axis [%g %g]',...
        mfilename,sprintf(' %g',sz),o.gr.zGr([end 1]));
end

if isempty(varargin)
    % ok dummy test
end

if isempty(sx) && isempty(sy) && isempty(sz)
    error('%s: you must give at least one of sx, sy or sz to define slices',mfilename);
end

% default axis properties
axProps = {'nextplot','add','xgrid','on','ygrid','on','zgrid','on'};
axProps = [axProps,'xlim',xlim,'ylim',ylim,'zlim',zlim];

%% Set up a nice figure as a template for the bubble environment
if ~isempty(figName)
    figure('name',figName,'position',figPos);
else
    if isempty(get(0,'children'))
        figure('name',figName,'position',figPos);
    else
        set(gcf,'position',figPos);
    end
end

wells{max(IComp)}  = NaN;
crange{max(IComp)} = NaN;
clim{max(IComp)}   = NaN;
ax(max(IComp))     = NaN;

k=0; small=1e-3;

for iComp= IComp(:)'
    if Delta
        well = well.setCout(deltaValues(o.UCN{iComp},STCONC{iComp}),iComp);
        crange{iComp} = ContourRange({o.UCN{iComp},STCONC{iComp}},o.numberOfConcContours);
        clim{iComp}   = crange{iComp}([1 end]);
    else
        well = well.setCout(o.UCN{iComp},iComp);
        crange{iComp} = ContourRange(o.UCN{iComp},o.numberOfConcContours);
        clim{iComp}   = crange{iComp}([1 end]);
    end

    k=k+1; % axes counter as IComop is not contiguous
    ax(iComp) = subplot(max(IComp),1,k,axProps{:},'clim',clim{iComp});
    view(3);
    
    % remove zero from crange (don't want a bubble of conc zero)
    crange{iComp} = crange{iComp}(abs(crange{iComp})>small);

    xlabel(ax(iComp),'x [m]');
    ylabel(ax(iComp),'y [m]');
    zlabel(ax(iComp),'z [m]');
    
        %% get googlemap
    try
        cprintf('Keywords','Trying to get a map from google, please wait ...');
        googleMap(ax(iComp),xlim,ylim,o.gr.zGr(end),map,coords);
        cprintf('Keywords',' got it.\n');
    catch ME
        cprintf(' failed.\n');
        warning(['mfLab:animate:googleMapFailed',...
                 '%s: %s\n',...
                 'Could not get google map ... I''ll continue without ...'],mfilename,ME.message);
    end

    % plot the grid at top
    if ~isempty(grdLineSpec)
        o.gr.plotGrid(ax(iComp),grdLineSpec);
    end

    if ~isempty(well)
        wells{iComp} = well.plot3D(ax(iComp));
    end

    light('Position',[-1 -1 1],'Style','infinite','parent',gca);
    lighting 'gouraud'
end

%% Simulation and video

time = o.t0+o.time;

set(gcf,'doublebuffer','on'); % guarantees flicker free video

[~,myname] = fileparts(mfilename); myname(1)=upper(myname(1));

vidObj           = VideoWriter([o.basename myname '.avi']);
vidObj.FrameRate = o.framerate;
vidObj.Quality   = o.quality;
vidObj.open();

%% Make video

ht(max(IComp)) = NaN;
pt=[];
for it=1:length(time)
    for iComp = IComp(:)'
        tts =[prefix o.titleStr{iComp} ', t = ', mf_datestr(time(it),o.dateformat,o.tdim)];
        
        if it==1
            ht(iComp) = title(ax(iComp),tts);
            if Delta
                hsl = slice(ax(iComp),o.gr.XM,o.gr.YM,o.gr.ZM,o.UCN{iComp}(it).values-STCONC{iComp},sx,sy,sz);
            else
                hsl = slice(ax(iComp),o.gr.XM,o.gr.YM,o.gr.ZM,o.UCN{iComp}(it).values,sx,sy,sz);
            end
            if ~isempty(well)
                wells{ iComp} = well.plot3D(ax(iComp));
            end
        else
            set(ht(iComp),'string',tts);
            if Delta
                hsl = slice(ax(iComp),o.gr.XM,o.gr.YM,o.gr.ZM,o.UCN{iComp}(it).values-STCONC{iComp},sx,sy,sz);
            else
                hsl = slice(ax(iComp),o.gr.XM,o.gr.YM,o.gr.ZM,o.UCN{iComp}(it).values,sx,sy,sz);               
            end
            if ~isempty(well)
                wells{ iComp} = wells{iComp}.plot3D(it);
            end
        end
    end
    
    set(ax(iComp),'clim',clim{iComp});
    
    set(hsl,'FaceColor','interp','EdgeColor','none',varargin{:},'facealpha',1);
    
    vidObj.writeVideo(getframe(gcf));
    
    if it>1
        delete(pt(pt>0));
    end
end

vidObj.close();
