function well = concXS(o,varargin)
% function animateObj = animateObj.concXS(gr,IComp,I,well,STCONC,varargin)
% plot XS either XS or YS of concentrations of the simiulation.
% gr = gridObj
% IComp = index/Nr of compoments/Species to plot
% I   = rows
% well= wells or MNW
% STCONC = STCONC may be omitted
% varargin = parameters to pass to the plotting routines and options
% options
%    'figure' figname
%    'figPos' figpos [xLL yLL w h] in pixels
%    'STCONC' STCONC
%    'xlim' xlim
%    'zlim'|'ylim' zlim (vertical dimension to show
%   'VK',VK to compute and visualize the hydraulic resistance of layers.
%
% if omitted --> defaults
%  figure = gcf
%  figPos = screenPos(0.75)
%  STCONC = [];
%  xLim   = gr.xc([1 end]);  same as gr.xGr([1 end])
%  yLim   = gr.zc([end  1]); same as gr.zGr([end 1])
%  I   = rows of wells if omitted and well exists depending on dir.
%  I   = 1 if model has only one row
%  IComp = all components 1:NCOMP
%
% TO 120908 121017 130403


[gr  ,varargin]     = getNext(varargin,'gridObj',[]);
if isempty(gr)
    error('%s: gridObj missing, is required',mfilename);
end

[IComp  ,varargin] = getNext(varargin,'double',[]);
[IComp  ,varargin] = getNext(varargin,'logical',IComp);
if isempty(IComp)
    IComp=1:o.NCOMP;
elseif islogical(IComp)
    IComp = find(IComp);
end
    
[IRow   ,varargin] = getNext(varargin,'double',[]);
[IRow   ,varargin] = getNext(varargin,'logical',IRow);
if  ~isempty(IRow) && ~isvector(IRow)    
    error(['%s: you must specify the rows to contour as third argument, like this:\n',...
           'animate.%s(gr,IComp,IRow,well,STCONC,varargin)'],mfilename,mfilename);
    IRow=1:gr.Ny;
end

[well,varargin]     = getNext(varargin,'wellObj',[]);
[well,varargin]     = getNext(varargin,'MNW1Obj',well);
[well,varargin]     = getNext(varargin,'MNW2Obj',well);

[STCONC ,varargin] = getNext(varargin,'double',[]);
[STCONC ,varargin] = getNext(varargin,'cell',STCONC);

if isempty(IRow)
    if ~isempty(well)
        IRow = unique([well.iy]);
    else
        error(['%s: without wells, you must specify the rows to contour as third argument, like this:\n',...
           'animate.%s(gr,IComp,IRow,well,STCONC,varargin)'],mfilename,mfilename);
    end
end

Delta = ~isempty(STCONC);
if Delta && ~iscell(STCONC)
    prefix='Delta';
    STCONC = {STCONC};
else
    prefix='';
end

[figName,varargin] = getProp(varargin,'fig',[]);
[figPos ,varargin] = getProp(varargin,'pos',screenPos(0.75));

[xlim   ,varargin] = getProp(varargin,'xlim',gr.xm([1 end]));
[ylim   ,varargin] = getProp(varargin,'ylim',gr.zm([end,1]));

[xLbl   ,varargin] = getProp(varargin,'xlabel','x [m]');
[yLbl   ,varargin] = getProp(varargin,'ylabel','z [m]');

% Undocumented option
[VK, ~ ] = getProp(varargin,'double',[]);

if ~isempty(VK)
    hydraulicC = gr.DZ./VK;
end

%% new figure explicitly requested ? (use 'fig','figname' in varargin)
if ~isempty(figName)
    figure('name',figName,'position',figPos);
else
    if isempty(get(0,'children'))
        figure('name',figName,'position',figPos);
    else
        set(gcf,'position',figPos);
    end
end

% Verify STCONC
if Delta
    if numel(o.UCN) ~= numel(STCONC),
        error('%s: number of compoments in STCONC (%s) does not match that in animateObj (%d)',...
            mfilename,numel(STCONC),numel(o.UCN));
    end
end

%% set cOut of wells of they exist

axProps = {'nextplot','add','xlim',xlim,'ylim',ylim,};

crange{max(IComp)} = NaN;
clim{max(IComp)}   = NaN;
ax(max(IComp))     = NaN;
tts{max(IComp)}    = NaN;
cmap = colormap();

k=0;
for iComp = IComp(:)'
    
    tts{iComp} = ['conc of ' prefix  o.titleStr{iComp}];

    if Delta
        if ~all(size(STCONC{iComp}) == size(o.UCN{iComp}(1).values))
           error('%s: size of arrays in STCONC [%d %d %d] does not match that of o.UNC [%d %d %d]  in animate',...
               mfilename,size(STCONC{iComp}),size(o.UCN{iComp}));
        end    

        crange{iComp} = ContourRange({o.UCN{iComp},STCONC{iComp}},o.numberOfConcCountours);            
        clim{iComp}   = crange{iComp}([1 end]);
        if ~isempty(well)
            well = well.setCout(deltaValues(o.UCN{iComp},STCONC{iComp}),iComp);
        end    
    else
        if ~all(gr.size==size(o.UCN{iComp}(1).values))
            error('%s: size of grid [%d %d %d] does not match that of conc{1} [%d %d %d] in animateObj',...
            mfilename,gr.size,size(o.UCN{iComp}.values));
        end

        crange{iComp} = ContourRange(o.UCN{iComp},o.numberOfConcContours);            
        clim{iComp}   = crange{iComp}([1 end]);
        if ~isempty(well)
            well = well.setCout(o.UCN{iComp},iComp);
        end    
    end

    k=k+1;    
    ax(iComp) = subplot(max(IComp),1,k,axProps{:},'clim',clim{iComp},'color',cmap(1,:));
    
    xlabel(ax(iComp),xLbl);
    ylabel(ax(iComp),yLbl);
    
end

%% Set up movie
vidObj = VideoWriter([o.basename '.avi']);
vidObj.FrameRate = o.framerate;
vidObj.Quality   = o.quality;
vidObj.open();

%% make movie
ht = NaN(o.NCOMP,1);
hC = NaN(o.NCOMP,1);
hP = NaN(o.NCOMP,1);

time = o.t0+o.time;

wells{max(IComp)} = NaN;
dateformat = 'yyyy-mm-dd';

for it=1:length(time)
    
    for iComp=1:o.NCOMP
        ttts = [tts{iComp}  '. Average of rows [' sprintf(' %d',IRow)  ' ]'  sprintf(', t = %s',datestr(time(it),dateformat))];
        if it==1
            ht(iComp) = title(ax(iComp),ttts);

            % concentration
            if Delta
                [~,hC(iComp)] = contourf(ax(iComp),gr.xc,gr.zc,...
                    XS(mean(o.UCN{iComp}(it).values(IRow,:,:)-STCONC{iComp}(:,:,IRow),1)),crange{iComp},'edgecolor','none');
            else
                [~,hC(iComp)] = contourf(ax(iComp),gr.xc,gr.zc,...
                    XS(mean(o.UCN{iComp}(it).values(IRow,:,:),1)),crange{iComp},'edgecolor','none');
            end
            colorbar('peer',ax(iComp));

            % stream lines
            if ~isempty(o.B) && ~isempty(o.prange) && gr.Ny==1
                [~,hP(iComp)] = contour(ax(iComp),gr.xp,gr.zp,o.B(it).Psi,o.prange,'color',grey);
            end
            
            if gr.Ny==1 && exist('hydraulicC','var')
                for iL= 1:gr.Nlay
                    if mean(hydraulicC(1,:,iL))>50
                        patch([gr.xGr gr.xGr(end:-1:1)], ...
                              [squeeze(gr.ZGR(1,:,iL+1)) squeeze(gr.ZGR(1,end:-1:1,iL))],...
                               'w',...
                              'facecolor','w','edgecolor','w','facealpha',0.25,'parent',ax(iComp));
                    end
                end          
            end

            % wells
            % different wells because of different axes with their own handles
            if ~isempty(well)
                wells{iComp} = well.plotXS(ax(iComp),IRow);
            end

        else
            set(ht(iComp),'string',ttts);
            
            if Delta
                set(hC(iComp),'zdata' ,XS(mean(o.UCN{iComp}(it).values(IRow,:,:)-STCONC{iComp}(IRow,:,:),1)));
            else
                set(hC(iComp),'zdata' ,XS(mean(o.UCN{iComp}(it).values(IRow,:,:),1)));
            end
            
            % stream lines
            if ~isempty(o.prange) && gr.Ny==1 && ~isempty(o.B)
                set(hP(iComp),'zdata',o.B(it).Psi);
            end
            
            % wells
            if ~isempty(well)
                wells{iComp} = wells{iComp}.plotXS(it,IRow);
            end
        end
    end
            
    vidObj.writeVideo(getframe(gcf));

end

vidObj.close();
