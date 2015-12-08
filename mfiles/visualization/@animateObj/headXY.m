function o = headXY(o,varargin)
%ANIMATEOBJ/HEADXY -- plot heads in the horizontal plane
%
% USAGE:
%    animate  = animateObj.headXY([gr,],ILay[,well[,STRTHD[,varargin]]]);
%
%
%  gr= gridObj;
%  gr       = gridObj
%  Ilay     = layers that will be averaged before computing contours
%  well     = wells or MNW (multi node wells, of class MNW1Obj)
%  STCONC   = STCONC may be omitted (change of conc will be plotted if present)
%  varargin may contain property name,value pairs.
%
%  if options are omitted --> defaults will be used:
%
%  recognized property,value pairs
%     'fig',figname
%     'pos',[x0 y0 dx dy], {screenPos(0.75)} ]
%     'plotgrid' [ {'c'} | linespec ]
%     'xlim', [ {gr.xGr([1 end]) } |  [xMin xMax] ]
%     'ylim', [ {gr.yGr([end 1]) } |  [yMin yMax] ]
%     'map',  [ {roadmap} | salellite | hybrid    ]
%     'co...' [ {rd} | wgs ]   % Notice: rd is rijksdriehoek (Dutch), wgs=google
%     'Ilay', [ {1 } | [validLayerNumbers(s)] [logial array for layers] ]
%     'mask',[mask1 mask1]
%     'drange',[range of value] contour elevations
%
%
% EXAMPLES:
%   animate = animate.headXY(gr,well,'iLay',5,'xlabel','x [m]','ylim',[10
%   animate = animate.headXY(gr,'iLay',5,'xlabel','x [m]','ylim',[10 50],...
%                    'fig','myfig','pos',[1 1 400 400],'color','ybgk',....)
% you may also use an absolute z value, 'z',-100 to select the layer to be
% plotted. Finally, you may also use an n by 3 array to define RGB colors.
%
% TO 120908 121017 130403

% TO 121219

if isempty(o.H)
    error('%s: animateObj does not contain heads, see call to animateObj to read heads',mfilename);
end

[o.gr      ,varargin] = getType(varargin,'gridObj',o.gr);
[~       ,varargin]   = getProp(varargin,'STRTHD',[]); % ignore STRTHD
[figName ,varargin] = getProp(varargin,'fig',[]);
[figPos  ,varargin] = getProp(varargin,'pos',screenPos(0.75));
[plgr    ,varargin] = getProp(varargin,'plotgrid',[]);
[z       ,varargin] = getProp(varargin,'z',[]);

if ~isempty(z), Ilay = find(o.gr.zGr<z,1,'first'); end

% Google map type ?
[map     ,varargin] = getProp(varargin,'map','roadmap');

% Coordinate system ? % use 'coo*','rd' in varargin
[coords  ,varargin] = getProp(varargin,'coord','rd');

if isempty(o.gr)
    error('%s: gridObj missing, is required because animateObj does not yet have one',mfilename);
end

[Ilay,varargin]     = getNext(varargin,'double',[]);
[Ilay,varargin]     = getNext(varargin,'logical',Ilay);
if isempty(Ilay), Ilay=1:o.gr.Nlay; end

[well,varargin]     = getNext(varargin,'wellObj',[]);
[well,varargin]     = getNext(varargin,'MNW1Obj',well);
[well,varargin]     = getNext(varargin,'MNW2Obj',well);

[~   ,varargin]     = getNext(varargin,'double',[]);  % ignore STRTHD


[xlim   ,varargin]  = getProp(varargin,'xlim',o.gr.xm([1 end]));
[ylim   ,varargin]  = getProp(varargin,'ylim',o.gr.ym([end,1]));

try
    time = [o.H.totim]+o.t0;
catch ME
    error(['%s: %s\n',...
        'Can''t find heads H in animateObj, call animateObj with {... ,''head'', ...}\n'],...
        mfilename,ME.message);    
end

% default axis properties
axProps = {'nextplot','add','xgrid','on','ygrid','on','zgrid','on','xlim',xlim,'ylim',ylim};

%% Set up a nice figure as a template for the bubble environment
if ~isempty(figName)
    figure('name',figName,'position',figPos);
else
    if isempty(get(0,'children'))
        figure('position',figPos);
    else
        set(gcf,'position',figPos);
    end
end

ax = axes(axProps{:});

%% Try to get a googlemap figure
try
   % googleMap, put it at bottom plane of drawing box
    googleMap(ax,xlim,ylim,o.gr.zGr(end),map,coords,'axis',ax,'rgb','on');
catch ME
    warning('on','mfLab:animate:googleMapFailed');
    warning('mfLab:animate:googleMapFailed',...
        '%s: %s\nCould not get google map ... I''ll continue without ...',mfilename,ME.message);
    warning('off','mfLab:animate:googleMapFailed');
end

%% Plot the grid ??
if ~isempty(plgr)
    if isLineSpec(plgr)
        o.gr.plotGrid(plgr);
    else
        o.gr.plotGrid('color','c','linewidth',0.1);
    end
end

[xLbl,varargin] = getProp(varargin,'xlabel','x [m]');
[yLbl,varargin] = getProp(varargin,'ylabel','y [m]');
[zLbl,varargin] = getProp(varargin,'zlabel','z [m]');

xlabel(ax,xLbl);
ylabel(ax,yLbl);
zlabel(ax,zLbl);

axis equal;
axis tight;

%% Set up video

[~,myname] = fileparts(mfilename); myname(1) = upper(myname(1));

set(gcf,'Renderer','OpenGL')

vidObj           = VideoWriter([o.basename myname]);
vidObj.FrameRate = o.framerate;
vidObj.Quality   = o.quality;
vidObj.open();

%% Make video
tt1 = 'Head';
[range  ,varargin] = getProp(varargin,'drange',[]);
if isempty(range)
    [range  , varargin] = getProp(varargin,'hrange',o.hrange);
end

[wellClr, varargin] = getProp(varargin,'color', 'k');

set(gca,'clim',range([1 end]));


for it=1:length(time)
        tts = sprintf('%0s range = [%.2f:%.2f:%.2f], t = %s',...
            tt1,range(1),range(end)-range(end-1),range(end),...
            datestr(time(it),'yyyy-mmm-dd')); 

        if it==1
            ht     = title(tts);
            [~,hc] = contour(ax,o.gr.xc,o.gr.yc,mean(o.H(it).values(:,:,Ilay),3),range,varargin{:});
            
%            set(get(hc,'children'),'CDataMapping','Direct');
            if ~isempty(well)
                well = well.plotXY(ax,'color',wellClr);
                well = well.plot3D(ax,'color',wellClr);                
            end
            colorbar('peer',ax);
        else
            set(ht,'string',tts);
            set(hc,'zdata',mean(o.H(it).values(:,:,Ilay),3));
            if ~isempty(well)
                well = well.plotXY(it);
            end
        end
            
    vidObj.writeVideo(getframe(gcf));    
end

vidObj.close();

fprintf('done\n');
