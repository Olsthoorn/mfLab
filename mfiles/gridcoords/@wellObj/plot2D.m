function o=plot2D(o,varargin)
    % well=well.plot2D(dir,ax|it[,I [,varargin]]]) -- plot or refresh plotted wells in cross section
    % well=well.plot2D(it);               -- update wells on 2D plot
    %
    % This function is best called from wellObj.plotXS and wellObj.plotYS.
    %
    % First call of this funcion must contain ax|it must be an axes.
    % subsequent calls contain ax|it should be the stress period number, it.
    %
    % I are the rows if dir=='x' or the columns if dir=='y'. The meaning is
    % that only those wells will be drawn that are in these rows or columns.
    %
    % dir = {'x'|'y'}
    % EXAMPLES
    %    well= well.plot2D('x',ax,Irow,varargin)  % plot on row iy
    %    well= well.plot2D('x',it,Irow,varargin)  % plot on row iy
    %    well= well.plot2D('y',it,Icol,varargin)  % plot on column ix
    %    well= well.plot2D(iper);               % refresh wells with current
    %    well= well.plot2D(ax,'y');
    % TO 120510 130330
    
    %% Perhaps we only need to refresh (assuming argument = it and whdl is a handle)
    
    if nargin<2
        error('%s: wellObj/plot2D must be called with axis or time index as first argument',mfilename);
    end
    
    if isempty(o)
        return;
    end
    
    % colors and alpha specified
    [EdgeColor ,varargin] = getProp(varargin,'edgeColor','');
    [FaceColor ,varargin] = getProp(varargin,'faceColor','');
    [FaceAlpha ,varargin] = getProp(varargin,'faceAlpha','');
    [EdgeAlpha ,varargin] = getProp(varargin,'edgeAlpha','');

    % if so, use it in all wells
    if ~isempty(EdgeColor), [o.EdgeColor] = deal(EdgeColor); end
    if ~isempty(FaceColor), [o.FaceColor] = deal(FaceColor); end
    if ~isempty(FaceAlpha), [o.FaceAlpha] = deal(FaceAlpha); end
    if ~isempty(EdgeAlpha), [o.EdgeAlpha] = deal(EdgeAlpha); end
    
    [dir,varargin]= getNext(varargin,'char','x');
    [ax,varargin] = getNext(varargin,'axis',[]);
    if ~isempty(ax)
        mustplot    = true;
        mustrefresh = false;
        it = 1;
    else
        mustplot    = false;
        mustrefresh = true;
        [it,varargin] = getNext(varargin,'double',1);        
    end

    [I,varargin ] = getNext(varargin,'double',[]);
    
    if ~isempty(I)
        for i=I(:)'
            for iw=1:numel(o)
                if dir=='x'
                    o(iw).UserData.mustPlot = ismember(i,o(iw).iy);
                else
                    o(iw).UserData.mustPlot = ismember(i,o(iw).ix);
                end
            end
        end
    else
        for iw=1:numel(o)
            o(iw).UserData.mustPlot = true;
        end
    end
        
% we dealt with handle, using variable "it" for refreshing
% item we may encouter in the varargin argument list
% error('wellObj/plot2D: property/value pairs must come in pairs.');
if mustrefresh
    for iw=1:length(o)
        if o(iw).UserData.mustPlot
            if isempty(o(iw).Q) ||it<1 || it>length(o(iw).Q)
            warning(['wewlObj:plot2D:outsideRange',...
                'wellObj/plot2D: Can''t refresh well because it=%d outside range[1:%d]'],...
                it,length(o(iw).Q));
                continue;
            end

            if ~all(ishandle(o(iw).whdl(1:2)))
                warning('wellObj:plot2D:propertiesNotinPairs',...
                    'wellObj/setWell: obj.whdl(1) must be a handle at this point, refreshing ignored');
                continue;
            end

            if isempty(o(iw).whdl) || ~strcmpi(get(o(iw).whdl(1),'type'),'patch')
                warning('wellObj:plot2D:wrongClassType',...
                    '%s: wellObj/setWell well.nr=%d: obj.whdl must be of type patch, may be empty.',mfilename,o(iw).nr);
                continue;
            end
            % now it is safe to refresh
            if  isnan(o(iw).Q(it)),
                set(o(iw).whdl(1:2),'visible','off');
            elseif o(iw).Q(it)==0
                set(o(iw).whdl(1),'visible','on');
                set(o(iw).whdl(2),'visible','off');
            else
                set(o(iw).whdl(1:2),'visible','on');                    
                set(o(iw).whdl(2),'facecolor',o(iw).FaceColor);
            end
    %         fprintf('well(Nr=%d), x=%.0f Q=%.0f whdl = [%s %s]\n',...
    %             o(iw).nr,o(iw).x,o(iw).Q(it),get(o(iw).whdl(1),'visible'),get(o(iw).whdl(2),'visible'));
        end
    end
    
    return;
end

if mustplot % only the selected wells given by logical array IW
    
    for iw= numel(o):-1:1
        if o(iw).UserData.mustPlot
            %% Plotting on this axis
            pxsz=get(gcf,'position');  % pixelsize
            if dir=='x'
                axWidth = diff(get(ax,'xlim')); % scale of the xaxis
                p       = o(iw).x;
            else
                axWidth = diff(get(ax,'xlim')); % not ylim
                p       = o(iw).y;
            end

            d=o(iw).wpix*axWidth/pxsz(3);      % world size of one pixel

            % d is the world width of 1 pixel, Wells will be drawn 2 pixels
            % wide minimum or as wide as twice the well radius if that
            % is more than the size of two pixels.
            if ~isempty(o(iw).rw), d=max(d,o(iw).rw); end

%             if length(varargin)>=1,
%                 o(iw).FaceColor=varargin{1}; %#ok<*AGROW>
%                 varargin(1)= [];
%             else
%                 o(iw).FaceColor='w';
%             end

%             if length(varargin)==2,
%                 o(iw).EdgeColor=varargin{2};
%                 varargin(1:2)=[];
%             else
%                 o(iw).EdgeColor='k';
%             end

            if ~isnan(o(iw).Q(1)),
                extra = {'visible','on'};
            else
                extra = {'visible', 'off'};
            end

            p=p(:)';
            zt = o(iw).ztop;
            z  = o(iw).z(:)';

            % casing above the screen            
            o(iw).whdl(1) = fill([p([1 1])-d, p([1 1])+d],[zt z(1) z(1) zt],grey,...
                'EdgeColor',o(iw).EdgeColor,'EdgeAlpha',o(iw).EdgeAlpha','parent',ax,extra{:}); % casing

            % screen (all points)
            o(iw).whdl(2) = fill([p-d, p(end:-1:1)+d],[z z(end:-1:1)],o(iw).FaceColor,...
                'EdgeColor',o(iw).EdgeColor,'FaceAlpha',o(iw).FaceAlpha,'EdgeAlpha',o(iw).EdgeAlpha,...
                'parent',ax,extra{:}); % screen

            set(o(iw).whdl(1),'facecolor','none');

            % remaining options are dealt with further down
            % using what is left in varargin

            if ~isempty(varargin)
                if rem(length(varargin),2) ~=0
                    warning('wellObj:plot2D:propertiesNoInPairs',...
                        'wellObj/plot2D: plotproperties/values must come in pairs.');
                else
                    for i=1:2:length(varargin)
                        if ~ischar(varargin{i})
                            warning('wellObj:plot2D:properiesNotChar',...
                                'wellObj/plot2D: plotproperties must be character strings');
                            continue;
                        else
                            set(o(iw).whdl(1),varargin{:});
                            set(o(iw).whdl(2),varargin{:});
                        end
                    end
                end
            end
        end
    end
end

