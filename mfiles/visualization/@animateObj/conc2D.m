function o = conc2D(o,varargin)
%ANIMATEOBJ/CONC2D -- shows transient conc in cross section along dir axis,
%where dir is 'x' or 'y'
%
% USAGE:
%    animateObj = animateObj.conc2D([gr,],dir,IComp,I,well,STCONC,varargin)
%
% gr       = gridObj
% IComp    = index/Nr of compoments/Species to plot
% I        = rows (average conc of the selected rows will be used)
% well     = wells or MNW (multi node wells, of class MNW1Obj)
% STCONC   = STCONC may be omitted
% varargin = parameters to pass to the plotting routines and options
% options that can be used as option,value pairs:
%    'figure' figname
%    'figPos' figpos [xLL yLL w h] in pixels
%    'STCONC' STCONC
%    'xlim' xlim
%    'zlim' zlim (vertical dimension to show)
%    'VK',VK to compute and visualize the hydraulic resistance of layers.
%
% if options are omitted --> defaults will be used:
%  figure = gcf
%  figPos = screenPos(0.75), that is 75% of the screen size
%  STCONC = [];
%  xLim   = gr.xc([1 end]);  same as gr.xGr([1 end])
%  yLim   = gr.zc([end  1]); same as gr.zGr([end 1])
%  I   = rows of wells to be shown, if omitted and if well objects exist
%  I   = 1 if model has only one row
%  IComp = all components 1:NCOMP
%  dir = 'x'
%
% TO 120908 121017 130403

%[topFig ,varargin] = getProp(varargin,'topFig'        ,screenPos(0.75));
[figName,varargin] = getProp(varargin,{'fig','figure'},'');
[figPos ,varargin] = getProp(varargin,{'figPos','pos'},screenPos(0.75));
[VK     ,varargin] = getProp(varargin,'double'        ,[]); % undocumented option  
[orient ,varargin] = getProp(varargin,'orien','V');  % figure orientation
if orient~='V'
    if ischar(orient)
        orient = upper(orient(1));
    else
        error('%s: option orien (figure orientation) must be ''V'' or ''H''',mfilename);
    end
end

[xLbl   ,varargin] = getProp(varargin,'xlabel','x [m]');
[yLbl   ,varargin] = getProp(varargin,'ylabel','y [m]');
[zLbl   ,varargin] = getProp(varargin,'zlabel','z [m]');

[dir  ,varargin] = getNext(varargin,'char','x');

if ~isempty(o.gr); % ignore gr in input
    [~,varargin] = getNext(varargin,'gridObj',[]);
else   % then get the grid
    [o.gr,varargin] = getNext(varargin,'gridObj',[]);
    if isempty(o.gr)
        error('%s: gridObj missing, is required because animateObj does not yet have one',mfilename);
    end
end

[xlim   ,varargin] = getProp(varargin,'xlim',o.gr.xc([1 end]));
[ylim   ,varargin] = getProp(varargin,'ylim',o.gr.yc([end,1]));
[zlim   ,varargin] = getProp(varargin,'zlim',o.gr.zc([end,1]));

[comp   ,varargin] = getProp(varargin,{'species','component'},[]);
if ~isempty(comp)
    IComp = strmatchi(comp,o.titleStr);
    if ~IComp(1)
        error('%s: unknown species <<%s>>, know species are %s',...
            mfilename,comp,sprintfs(' %s',o.titleStr));
    end
else
    [comp   ,varargin] = getNext(varargin,{'cell','char'},[]);
    if ~isempty(comp)        
        IComp = strmatchi(comp,o.titleStr);
        if ~IComp(1)
            error('%s: animateObj does not have this/these compoments <<%s>>, it has <<%s>>',...
                mfilename,sprintfs(' %s',comp),sprintf(' %s',o.titleStr));
        end
    else
        [IComp  ,varargin] = getNext(varargin,'double',[]);
        [IComp  ,varargin] = getNext(varargin,'logical',IComp);
    end
end

if     isempty(IComp),     IComp = 1:o.NCOMP;
elseif islogical(IComp),   IComp = find(IComp);
end
    
[I   ,varargin] = getNext(varargin,'double',[]);
[I   ,varargin] = getNext(varargin,'logical',I);

[well,varargin]     = getNext(varargin,'wellObj',[]);
[well,varargin]     = getNext(varargin,'MNW1Obj',well);
[well,varargin]     = getNext(varargin,'MNW2Obj',well);


%% without wells you must specify I
if isempty(I)
    if dir=='x'
        if o.gr.Ny==1
            I = 1;
        else
            if ~isempty(well)
                I = [well.iy];
            else
                error(['%s: without wells, you must specify the rows to contour as third argument, like this:\n',...
                   'animate.%s(gr,IComp,I,well,STCONC,varargin)'],mfilename,mfilename);
            end
        end
    else
        if o.gr.Nx==1
            I = 1;
        else
            if ~isempty(well)
                I = [well.ix];
            else
                error(['%s: without wells, you must specify the columns to contour as third argument, like this:\n',...
                   'animate.%s(gr,IComp,I,well,STCONC,varargin)'],mfilename,mfilename);
            end
        end
    end
end

if  ~(isvector(IComp) && isvector(I)) || isempty(I)
    error(['%s: I probably read STCONC in place of IComp or IRowCol.\n',...
           'REMEDY: specify the IComp and IRowCol explicitly and in the right\n',...
           'order. See help %s to find out the right order.\n',...
           'animate.%s(gr,IComp,IRowCol,well,STCONC,varargin)'],mfilename,mfilename);
end

[STCONC ,varargin] = getNext(varargin,'cell',[]);
[STCONC ,   ~    ] = getNext(varargin,'double',STCONC);

Delta = ~isempty(STCONC);
if Delta && ~iscell(STCONC)
    prefix='Delta';
    STCONC = {STCONC};
else
    prefix='';
end

% Undocumented option

if ~isempty(VK)
    hydraulicC = gr.DZ./VK;
end

%% new figure explicitly requested ? (use 'fig','figname' in varargin)
if ~isempty(figName)
    figure('name',figName,'position',figPos);
else
    figure('name',pwd,'position',figPos);
end

% Verify STCONC
if Delta
    if numel(o.UCN) ~= numel(STCONC),
        error('%s: number of compoments in STCONC (%s) does not match that in animateObj (%d)',...
            mfilename,numel(STCONC),numel(o.UCN));
    end
end

%% set cOut of wells of they exist

crange{max(IComp)} = NaN;
clim{max(IComp)}   = NaN;
ax(max(IComp))     = NaN;
cmap = colormap();

k=0;
for iComp = IComp(:)'

    if Delta
        if ~all(size(STCONC{iComp}) == size(o.UCN{iComp}(iComp).values))
           error('%s: size of arrays in STCONC [%d %d %d] does not match that of o.UNC [%d %d %d]  in animate',...
               mfilename,size(STCONC{iComp}),size(o.UCN{iComp}));
        end    

        crange{iComp} = ContourRange({o.UCN{iComp},STCONC{iComp}},o.numberOfConcContours);            
        clim{iComp}   = crange{iComp}([1 end]);
        if ~isempty(well)
            well = well.setCout(deltaValues(o.UCN{iComp},STCONC{iComp}),iComp);
        end    
    else
        if ~all(o.gr.size==size(o.UCN{iComp}(1).values))
            error('%s: size of grid [%d %d %d] does not match that of conc{1} [%d %d %d] in animateObj',...
            mfilename,o.gr.size,size(o.UCN{iComp}.values));
        end

        crange{iComp} = ContourRange(o.UCN{iComp},o.numberOfConcContours);
        if isempty(crange{iComp})
            error('%s: Nothing to contour for compoment %d = ''%s''',...
                mfilename,iComp,o.titleStr{iComp});
        end
        clim{iComp}   = crange{iComp}([1 end]);
        if ~isempty(well)
            well = well.setCout(o.UCN{iComp},iComp);
        end    
    end

    k=k+1;
    if orient=='H'
        ax(iComp) = subplot(1,max(IComp),k,'clim',clim{iComp},'color',cmap(1,:));
    else
        ax(iComp) = subplot(max(IComp),1,k,'clim',clim{iComp},'color',cmap(1,:));
    end
    
end

set(ax,'fontsize',14);

%% Set up movie
vidObj = VideoWriter([o.basename '.avi']);
vidObj.FrameRate = o.framerate;
vidObj.Quality   = o.quality;
vidObj.open();

%% make movie
ht = NaN(o.NCOMP,1);
hC = NaN(o.NCOMP,1);
hP = NaN(o.NCOMP,1);

aYear = 365.24;

if o.t0>datenum(1000,1,1)
    time = o.t0+o.time;
    o.tdim='AD';
elseif o.time(end) > 3000
    time = o.t0+o.time;
    o.tdim='y';
else
    time = o.t0+o.time;
    o.tdim='d';
end


wells{max(IComp)} = NaN;

for it=1:length(time)
    
    switch o.tdim
        case 'd'
            timeStr = sprintf('%.0f d',time(it));
        case 'AD'
            timeStr = datestr(time(it),o.dateformat);
        case 'y'
            timeStr = sprintf('%.0f y',time(it)/aYear);
        otherwise
    end
    
    for iComp=IComp(:)'
        tts = [prefix o.titleStr{iComp}  '. Average of rows [' sprintf(' %d',I)  ' ], t = ' timeStr];
        if it==1
            ht(iComp) = title(ax(iComp),tts,'fontsize',14);

            % concentration
            if dir=='x'
                set(ax(iComp),'nextplot','add','xlim',xlim,'ylim',zlim);

                if Delta
                    [~,hC(iComp)] = contourf(ax(iComp),o.gr.xc,o.gr.zc,...
                        XS(mean(o.UCN{iComp}(it).values(I,:,:)-STCONC{iComp}(I,:,:),1)),crange{iComp},'edgecolor','none');
                else
                    [~,hC(iComp)] = contourf(ax(iComp),o.gr.xc,o.gr.zc,...
                        XS(mean(o.UCN{iComp}(it).values(I,:,:),1)),crange{iComp},'edgecolor','none');
                end
                
                % stream lines
                if ~isempty(o.B) && ~isempty(o.prange) && o.gr.Ny==1
                    [~,hP(iComp)] = contour(ax(iComp),o.gr.xp,o.gr.zp,o.B(it).Psi,o.prange,'color',grey);
                end
                
                if o.gr.Ny==1 && exist('hydraulicC','var')
                    for iL= 1:o.gr.Nlay
                        if mean(hydraulicC(1,:,iL))>50
                            patch([o.gr.xGr o.gr.xGr(end:-1:1)], ...
                                  [squeeze(o.gr.ZGR(1,:,iL+1)) squeeze(o.gr.ZGR(1,end:-1:1,iL))],...
                                   'w',...
                                  'facecolor','w','edgecolor','w','facealpha',0.25,'parent',ax(iComp));
                        end
                    end          
                end

                % wells
                % different wells because of different axes with their own handles
                if ~isempty(well)
                    wells{iComp} = well.plotXS(ax(iComp),I);
                end

                xlabel(ax(iComp),xLbl,'fontsize',14);
                ylabel(ax(iComp),zLbl,'fontsize',14);
                                
            else
                set(ax(iComp),'nextplot','add','xlim',ylim,'ylim',zlim);

                if Delta
                    [~,hC(iComp)] = contourf(ax(iComp),o.gr.yc,o.gr.zc,...
                        YS(mean(o.UCN{iComp}(it).values(:,I,:)-STCONC{iComp}(:,I,:),2)),crange{iComp},'edgecolor','none');
                else
                    [~,hC(iComp)] = contourf(ax(iComp),o.gr.yc,o.gr.zc,...
                        YS(mean(o.UCN{iComp}(it).values(:,I,:),2)),crange{iComp},'edgecolor','none');
                end
                
                % stream lines
                if ~isempty(o.B) && ~isempty(o.prange) && o.gr.Nx==1
                    if ~isfield(o.B,'PhiY')
                        o.B = mf_Psi(o.B,1,'y');
                    end
                    [~,hP(iComp)] = contour(ax(iComp),o.gr.yp,o.gr.zp,o.B(it).PsiY,o.prange,'color',grey);
                end

                if o.gr.Nx==1 && exist('hydraulicC','var')
                    for iL= 1:o.gr.Nlay
                        if mean(hydraulicC(1,:,iL))>50
                            patch([o.gr.yGr o.gr.yGr(end:-1:1)], ...
                                  [squeeze(o.gr.ZGR(1,:,iL+1)) squeeze(o.gr.ZGR(1,end:-1:1,iL))],...
                                   'w',...
                                  'facecolor','w','edgecolor','w','facealpha',0.25,'parent',ax(iComp));
                        end
                    end          
                end

                % wells
                % different wells because of different axes with their own handles
                if ~isempty(well)
                    wells{iComp} = well.plotYS(ax(iComp),I);
                end

                xlabel(ax(iComp),yLbl);
                ylabel(ax(iComp),zLbl);
            end
            
            if ~isfield(o.patch,'options')
                for iP=1:numel(o.patch)                    
                    o.patch(iP).options = {};
                end
            end

            if  ~isempty(o.patch)
                for iP=1:numel(o.patch)
                    patch(o.patch(iP).x, o.patch(iP).y,o.patch(iP).color,...
                        o.patch(iP).options{:});
                end
            end
            
% next line does not work in Matlab 2014b                
%            set(get(colorbar('peer',ax(iComp)),'title'),'string',o.titleStr{iComp}); 
% so we use these 3 steps instead
            axes(ax(iComp)); %#ok switch focus to ax(iComp)
            cAx = colorbar;  % use current fig (no arguments)
            set(get(cAx,'title'),'string',o.titleStr{iComp}); % then set title

        else
            set(ht(iComp),'string',tts);
            
            if dir=='x'
                if Delta
                    set(hC(iComp),'zdata' ,XS(mean(o.UCN{iComp}(it).values(I,:,:)-STCONC{iComp}(I,:,:),1)));
                else
                    set(hC(iComp),'zdata' ,XS(mean(o.UCN{iComp}(it).values(I,:,:),1)));
                end
                
                % stream lines
                if ~isempty(o.prange) && o.gr.Ny==1 && ~isempty(o.B)
                    set(hP(iComp),'zdata',o.B(it).Psi);
                end

                % wells
                if ~isempty(well)
                    wells{iComp} = wells{iComp}.plotXS(it,I);
                end

            else
                if Delta
                    set(hC(iComp),'zdata' ,YS(mean(o.UCN{iComp}(it).values(:,I,:)-STCONC{iComp}(:,I,:),2)));
                else
                    set(hC(iComp),'zdata' ,YS(mean(o.UCN{iComp}(it).values(:,I,:),2)));
                end
                
                % stream lines
                if ~isempty(o.prange) && o.gr.Ny==1 && ~isempty(o.B)
                    set(hP(iComp),'zdata',o.B(it).PsiY);
                end
                
                % wells
                if ~isempty(well)
                    wells{iComp} = wells{iComp}.plotYS(it,I);
                end

            end            
            
        end
    end
            
    vidObj.writeVideo(getframe(gcf));

end

vidObj.close();
