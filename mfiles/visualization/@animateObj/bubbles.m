function o = bubbles(o,varargin)
 % aninmateObj.bubles([gr,] IComp,well,STCONC)
 % animate bubbles if STCONC is given, subtract them first
 %
 % TO 120914

[figPos ,varargin] = getProp(varargin,{'figpos','pos}',screenPos(0.75));
[figName,varargin] = getProp(varargin,'fig','');
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
if isempty(well)
    error('animateObj:bubles:noWellsPresent',...
        ['%s: method %s requires well or multinode wells in call.\n',...
          'REMEDY: call like     animate = animateObj.%s(gr,well,...)'],...
          mfilename,mfilename,mfilename);
end

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
[zlim, ~       ] = getProp(varargin,'zlim',o.gr.zGr([end 1]));

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

%% Set a convenient maximum limit of bubble contours
numberOfBubbleContours = 10;
bubbleClr  = 'br';  % blue and red

wells{max(IComp)}  = NaN;
crange{max(IComp)} = NaN;
clim{max(IComp)}   = NaN;
ax(max(IComp))     = NaN;

k=0; small=1e-3;

for iComp= IComp(:)'
    if Delta
        well = well.setCout(deltaValues(o.UCN{iComp},STCONC{iComp}),iComp);
        crange{iComp} = ContourRange({o.UCN{iComp},STCONC{iComp}},numberOfBubbleContours);
        crange{iComp} = -2.5:1:2.5;
        clim{iComp}   = crange{iComp}([1 end]);
    else
        well = well.setCout(o.UCN{iComp},iComp);
        crange{iComp} = ContourRange({o.UCN{iComp},STCONC{iComp}},numberOfBubbleContours);
        clim{iComp}  =  -2.5:1:2.5;
        clim{iComp}   = crange{iComp}([1 end]);
    end

    k=k+1; % axes counter as IComop is not contiguous
    ax(iComp) = subplot(max(IComp),1,k,axProps{:});
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

    % plot teh grid at top
    if ~isempty(grdLineSpec)
        o.gr.plotGrid(ax(iComp),grdLineSpec);
    end

    if ~isempty(well)
        wells{iComp} = well.plot3D(ax(iComp));
    end

    light('Position',[-1 -1 1],'Style','infinite','parent',gca);
    lighting 'gouraud'
end

% %% We will use as many axes on the figure as we have components.
% % We do this by copying our marked up teplate (ax on hfig)
% 
% % Parameters (margins) to compute position of the axes:
% oLB =0.1; oM=0.08; oRT=0; % leftBOttom, mid and rightTop relative margins
% ih=1; nh = round(sqrt(o.NCOMP)); dh=(1-oLB-nh*oM-oRT)/nh;  % hor position of axis
% iv=1; nv = ceil( sqrt(o.NCOMP)); dv=(1-oLB-nv*oM-oRT)/nv;  % ver position of axis
% 
% % Copy and repostion axis for all NCOMP species/components
% ax = NaN(1,o.NCOMP); % allocate axes
% 
% for iComp= 1:o.NCOMP
%     axpos = [oLB+(ih-1)*(dh+oM), oLB+(iv-1)*(dv+oM), dh, dv];
%     ih = ih+1; if ih>nh, ih=1; iv=iv+1; end
%     
%     if iComp==1,
%         ax(iComp) = ax;                
%     else
%         ax(iComp) = copyobj(ax(1),get(ax(1),'parent'));
%     end
%     set(ax(iComp),'position',axpos);
% end

%% Simulation and video

time = o.t0+o.time;

%% Set up video
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
    for iComp = 1:o.NCOMP
        tts =[prefix o.titleStr{iComp} ', isosurf [' sprintf(' %g',crange{iComp}) '] t = ', mf_datestr(time(it),o.dateformat.o.tdim)];
        
        if it==1
            ht(iComp) = title(ax(iComp),tts);
            %wells{ iComp} = wells{iComp}.plot3D(ax(iComp));
        else
            set(ht(iComp),'string',tts);  
            % wells must exist for bubbles, had be asserted before
            wells{ iComp} = wells{iComp}.plot3D(it);
        end
            
        for iSurf = numel(crange{iComp}):-1:1
            if Delta
                C = o.UCN{iComp}(it).values - STCONC{iComp};
                clr = bubbleClr((crange{iComp}(iSurf)>0) + 1);
            else
                C = o.UCN{iComop}(it).values;
                clr = bubbleClr((crange{iComp}(iSurf)>mean(crange{iComp})) + 1);
            end
            
            FV = isosurface(o.gr.XM,o.gr.YM,o.gr.ZM,C,crange{iComp}(iSurf));
            pt(iComp,iSurf) = patch(FV,'parent',ax(iComp),...
                'facealpha',0.25,'facecolor',clr,'edgecolor','none');
            isonormals(o.gr.XM,o.gr.YM,o.gr.ZM,C,pt(iComp,iSurf));
            
        end                
    end
    vidObj.writeVideo(getframe(gcf));
    
    if it>1
        delete(pt(pt>0));
    end
end

vidObj.close();
