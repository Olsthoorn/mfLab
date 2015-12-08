function o = head2D(o,varargin)
%ANIMATEOBJ/HEAD2D -- shows transient heads in cross section along dir-axis
%
% USAGE:
%    animateObj = animateObj.head2D(gr,dir,I,well,STRTHD,varargin)
%
% gr       = gridObj
% IComp    = index/Nr of compoments/Species to plot
% I        = row number if dir='x' or column number of dir='y'  or 1 if omitted
% well     = wells or MNW (multi node wells, of class MNW1Obj)
% STRTHD   = STRTHD may be omitted (causes drawdown to be plotted if present)
% varargin = parameters to pass to the plotting routines and options
% options that can be used as option,value pairs:
%    'figure' figname
%    'figPos' figpos [xLL yLL w h] in pixels
%    'STRTHD' STRTHD
%    'xlim' xlim
%    'zlim' zlim (vertical dimension to show)
%    'VK',VK to compute and visualize the hydraulic resistance of layers.
%
% if options are omitted --> defaults will be used:
%  figure = gcf
%  figPos = screenPos(0.75), that is 75% of the screen size
%  STRTHD = [];
%  xLim   = gr.xc([1 end]);  same as gr.xGr([1 end])
%  yLim   = gr.zc([end  1]); same as gr.zGr([end 1])
%  I   = rows if dir='x' or columns if dir='y' of wells to be shown, if omitted and if well objects exist
%  I   = 1 if model has only one row
%
% TO 120908 121017 130403

[figPos ,varargin] = getProp(varargin,{'figPos','pos'},screenPos(0.75));
[figName,varargin] = getProp(varargin,{'fig','figure'},'');

[o.gr,varargin] = getThis(varargin,'gridObj',o.gr);

[xlim   ,varargin] = getProp(varargin,'xlim',o.gr.xc([1 end]));
[ylim   ,varargin] = getProp(varargin,'ylim',o.gr.yc([end,1]));
[zlim   ,varargin] = getProp(varargin,'zlim',o.gr.zc([end,1]));


[xLbl   ,varargin] = getProp(varargin,'xlabel','x [m]');
[yLbl   ,varargin] = getProp(varargin,'ylabel','y [m]');
[zLbl   ,varargin] = getProp(varargin,'zlabel','z [m]');


[dir,varargin]  = getNext(varargin,'char','x');


[I  ,varargin]    = getNext(varargin,'double',[]);

[well,varargin] = getNext(varargin,'wellObj',[]);
[well,varargin] = getNext(varargin,'MNW1Obj',well);
[well,varargin] = getNext(varargin,'MNW2Obj',well);

if numel(size(I))>2
    error(['%s: Probably STRTHD was read in place of IRowCol.\n',...
        'IRowCol may be empty or a numerical vector, one of them,\n',...
        'was a matrix or 3D array. You have to check that and make sure that\n',...
        'the variables are in the right order. If necessary specify IRowCol explicitly'],...
        mfilename);
end

if isempty(I)
    if isempty(well)
        if dir == 'x'
                I = 1:o.gr.Ny;
        else
                I = 1:o.gr.Nx;
        end
    else
        if dir == 'x'
            I = unique([well.iy]);
        else
            I = unique([well.ix]);
        end
    end
else
    if ~isvector(I)
        error(['%s: I must be a vector, it is an array, you may have forgotton\n',...
            'to specify it explicitly and now I got STRTHD instead.',....
            'REMEDY: specify I explicitly to make the input unique'],...
            mfilename);
    end
end

%% Set all heads lower than their cell bottom equal to NaN
for it = 1:numel(o.H)
    o.H(it).values(o.H(it).values<o.gr.ZBlay)=NaN;
end

[STRTHD,varargin] = getNext(varargin,'double',[]);
[STRTHD,varargin] = getNext(varargin,'cell',STRTHD);

Delta = ~isempty(STRTHD);
if Delta && ~iscell(STRTHD)
    prefix='Delta';
    STRTHD = {STRTHD};
else
    prefix='';
end

%varargin(cellfun(@isempty,varargin))=[];

if o.gr.AXIAL, xLbl = 'r [m]'; end

%% New figure
if ~isempty(figName)
    figure('name',figName,'position',figPos);
else
    if isempty(get(0,'children'))
        figure('name',figName,'position',figPos);
    else
        set(gcf,'position',figPos);
    end
end

if dir =='x'
    ax = axes('nextplot','add','xlim',xlim,'ylim',zlim,'clim',o.hrange([1 end]));
    xlabel(ax,xLbl);
    ylabel(ax,zLbl);
else
    ax = axes('nextplot','add','xlim',ylim,'ylim',zlim,'clim',o.hrange([1 end]));
    xlabel(ax,yLbl);
    ylabel(ax,zLbl);
end

%% Set up movie
vidObj = VideoWriter([o.basename '.avi']);
vidObj.FrameRate = o.framerate;
vidObj.Quality   = o.quality;
vidObj.open();

%% make movie
o.time     = o.t0 + o.time;

for it=1:length(o.time)
    
    tts = [prefix 'Head, I = [' sprintf(' %d',I) ' ]' sprintf(', time = %s',mf_datestr(o.time(it),o.dateformat,o.tdim))];

    if it==1
        ht = title(ax,tts);

        if dir == 'x'
            if Delta
                [~,hH] = contourf(ax,o.gr.xc,o.gr.zc,...
                    XS(mean(o.H(it).values(I,:,:)-STRTHD(I,:,:),1)),o.hrange,...
                    'edgecolor','none',varargin{:});
            else
                [~,hH] = contourf(ax,o.gr.xc,o.gr.zc,...
                    XS(mean(o.H(it).values(I,:,:),1)),o.hrange,'edgecolor','none',varargin{:});
            end
            
            % stream lines
            if o.gr.Ny==1 && ~isempty(o.B)
                [~,hP] = contour(ax,o.gr.xp,o.gr.zp,o.B(it).Psi,o.prange,'color',grey,varargin{:});
            end

            % wells
            % different wells because of different axes with their own handles
            if ~isempty(well)
                well = well.plotXS(ax,I,'w');
            end

        else
            if Delta
                [~,hH] = contourf(ax,o.gr.yc,o.gr.zc,...
                    YS(mean(o.H(it).values(:,I,:)-STRTHD(:,I,:),2)),o.hrange,...
                    'edgecolor','none',varargin{:});
            else
                [~,hH] = contourf(ax,o.gr.yc,o.gr.zc,...
                    YS(mean(o.H(it).values(:,I,:),2)),o.hrange,'edgecolor','none',varargin{:});
            end
            
            % stream lines
            o.B = mf_Psi(o.B,I,'y');
            if o.gr.Nx==1 && ~isempty(o.B)
                [~,hP] = contour(ax,o.gr.yp,o.gr.zp,o.B(it).PsiY,o.prange,'color',grey,varargin{:});
            end

            % wells
            % different wells because of different axes with their own handles
            if ~isempty(well)
                well = well.plotYS(ax,I,'w');
            end
            
        end
        
        if  ~isempty(o.patch)
            for iP=1:numel(o.patch)
                patch(o.patch(iP).x, o.patch(iP).y,o.patch(iP).color,o.patch(iP).options{:});
            end
        end
        
        colorbar('peer',ax);
        
    else
        set(ht,'string',tts);

        if dir == 'x'
            if Delta
                set(hH,'zdata' ,XS(mean(o.H(it).values(I,:,:)-STRTHD(I,:,:),1)));
            else
                set(hH,'zdata' ,XS(mean(o.H(it).values(I,:,:),1)));
            end

            % stream lines
            if o.gr.Ny==1 && ~isempty(o.B)
                set(hP,'zdata',o.B(it).Psi);
            end

            % wells
            if ~isempty(well)
                well = well.plotXS(it);
            end
        else
            if Delta
                set(hH,'zdata' ,YS(mean(o.H(it).values(:,I,:)-STRTHD(:,I,:),2)));
            else
                set(hH,'zdata' ,YS(mean(o.H(it).values(:,I,:),2)));
            end

            % stream lines
            if o.gr.Nx==1 && ~isempty(o.B)
                set(hP,'zdata',o.B(it).PsiY);
            end

            % wells
            if ~isempty(well)
                well = well.plotYS(it);
            end
        end
    end
            
    vidObj.writeVideo(getframe(gcf));

end

vidObj.close();
