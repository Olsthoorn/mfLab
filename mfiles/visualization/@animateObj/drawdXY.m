function o = drawdXY(o,varargin)
% aninmateObj.headXY(gr[,well,STRTHD,varargin)
% animate drawd in plane XY
% example: using pvpairs (property value pairs)
% recognized property,value pairs
%     'fig', [{''} | figname }
%     'pos', [ {screenPos(0.75)} | [] ]
%     'plotgrid' [ {'c'} | linespec ]
%     'xlim', [ {gr.xGr([1 end]) } |  [xMin xMax] ]
%     'ylim', [ {gr.yGr([end 1]) } |  [yMin yMax] ]
%     'map',  [ {roadmap} | salellite | hybrid    ]
%     'co...' [ {rd} | wgs ]   % rd is rijksdriehoek (dutch), wgs=google
%     'iz',   [ {1 } | [validLayerNumbers(s)] ]
%     'z',    [ { gr.zLay(1) } | validZvalue ]
% animate = animate.headXY(gr[,well],'iz',5,'xlabel','x [m]','ylim',[10
% 50],'fig','myfig','pos',[1 1 400 400],'color','ybgk',....)
% you may also use an absolute z value, 'z',-100
% you may also use an [n,3] array of RGB values for color
%
% TO 121219

if isempty(o.H)
    error('%s: animageObj does not contain heads, see call to animateObj to read heads',mfilename);
end

[STRTHD,varargin] = getProp(varargin,'STRTHD',[]);
[figName ,varargin] = getProp(varargin,'fig',[]);
[figPos  ,varargin] = getProp(varargin,'pos',screenPos(0.75));
[plgr    ,varargin] = getProp(varargin,'plotgrid',[]);

% Google map type ?
[map     ,varargin] = getProp(varargin,'map','roadmap');

% Coordinate system ? % use 'coo*','rd' in varargin
[coords  ,varargin] = getProp(varargin,'coord','rd');

% Make sure varargin exists and never gets empty
[gr  ,varargin]     = getNext(varargin,'gridObj',[]);
[well,varargin]     = getNext(varargin,'wellObj',[]);
[well,varargin]     = getNext(varargin,'MNW1Obj',well);
[well,varargin]     = getNext(varargin,'MNW2Obj',well);
[STRTHD,varargin]   = getNext(varargin,'double',STRTHD);

if isempty(gr)
    error('%s: gridObj missing, is required',mfilename);
end

if ~isempty(o.D)
    % o.D will be used
else
    if isempty(o.H)
        error('%s: Neither heads nor drawdowns given, can''t contour drawdowns\n',mfilename);
    else
        if isempty(STRTHD)
            error(['%s: Need STRTHD as 3rd argument of %s to contour drawdowns\n',...
                'REMEDY: call animateObj with ''drawdown'' or call %s with STRTHD as 3rd argument'],...
                mfilename,mfilename,mfilename);
        else
            o.D = o.H;
            for it=1:numel(o.H)
                o.D(it).values = o.D(it).values-STRTHD;
            end
        end
    end
end

[xlim   ,varargin]  = getProp(varargin,'xlim',gr.xm([1 end]));
[ylim   ,varargin]  = getProp(varargin,'ylim',gr.ym([end,1]));

try
    time = [o.H.totim]+o.t0;
catch ME
    error(['%s: %s\n',...
        'Can''t find heads H in animateObj, call animateObj with {... ,''drawd'', ...}\n'],...
        mfilename,ME.message);    
end

% default axis properties
axprops = {'nextplot','add','xgrid','on','ygrid','on','zgrid','on','xlim',xlim,'ylim',ylim};

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

ax = axes(axprops{:});

%% Usem might have specified iz or z
[iz,varargin] = getProp(varargin,'iz',[]);
if isempty(iz)
    [z,varargin] = getProp(varargin,'z',[]);
    if isempty(z)
        iz = 1;
    else
        iz = hit(gr.zGr,z);
    end
end

%% Try to get a googlemap figure
try
   % googleMap, put it at bottom plane of drawing box
    googleMap(ax,xlim,ylim,gr.zGr(end),map,coords,'axis',ax,'rgb','on');
catch %#ok
    warning('mfLab:animate:googleMapFailed',...
        '%s: Could not get google map ... I''ll continue without ...',mfilename);
end

%% Plot the grid ??
if ~isempty(plgr)
    if isLineSpec(plgr)
        gr.plotGrid(plgr);
    else
        gr.plotGrid('color','c','linewidth',0.1);
    end
end

[xLbl,varargin] = getProp(varargin,'xlabel','x [m]');
[yLbl,varargin] = getProp(varargin,'ylabel','y [m]');
[zLbl,varargin] = getProp(varargin,'zlabel','z [m]');

xlabel(ax,xLbl);
ylabel(ax,yLbl);
zlabel(ax,zLbl);

axis equal;

%% Set up video

[~,myname] = fileparts(mfilename); myname(1) = upper(myname(1));

set(gcf,'Renderer','OpenGL')

vidObj           = VideoWriter([o.basename myname]);
vidObj.FrameRate = o.framerate;
vidObj.Quality   = o.quality;
vidObj.open();

%% Make video
tt1 = 'Drawdown';

[range  , varargin] = getProp(varargin,'drange',o.drange);
[wellClr, varargin] = getProp(varargin,'color','k');

set(gca,'clim',range([1 end]));

for it=1:length(time)
        tts = sprintf('%s range = [%.2f:%.2f:%.2f], t = %s',...
            tt1,range(1),range(end)-range(end-1),range(end),...
            datestr(time(it),'yyyy-mmm-dd')); 

        if it==1
            ht = title(ax,tts);
                [~,hc] = contour(ax,gr.xc,gr.yc,o.D(it).values(:,:,iz),range,varargin{:});
%            set(get(hc,'children'),'CDataMapping','Direct');
            well = well.plotXY(ax,'color',wellClr);
            
        else
            set(ht,'string',tts);
            set(hc,'zdata',o.D(it).values(:,:,iz));
            well = well.plotXY(it);
        end
            
    vidObj.writeVideo(getframe(gcf));    
end

vidObj.close();

fprintf('done\n');
