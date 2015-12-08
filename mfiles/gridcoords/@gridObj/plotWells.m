function well=plotWells(~,varargin)
% OBSOLETE
%% hw = gr.plotWells(ax,wells)  -- plots   wellscreens in 3D
%% hw = gr.plotWells(it,wells,hw)  -- updates wellscreens in 3D
%  based on well(iw).Q(it) being NaN, zero or non-zero
% 
% TO 120531

grey     = [0.8 0.8 0.8];
defaultcolor     = 'b';
defaultlinewidth = 2;
mustPlot = false;

if ~exist('varargin','var') || isempty(varargin)
    error('oidObj:plotWells:noData',...
        '%s: argument wellObj needed',mfilename);
end

try
    ax=varargin{1};
    if ~strcmp(get(ax,'type'),'axes')
        eror();
    end
    mustPlot=true;
    it = 1;
catch ME %#ok
    if ~isnumeric(varargin{1});
        error('oid:plotWells:noUpdateIndex',...
        '%s: update index needed for wellObj',mfilename);
    else
        it = varargin{1};
    end
end

varargin(1)=[];

if isempty(varargin)
    error('oid:plotWells:noWellData',...
        'oidObj/%s: argument wellObj needed',mfilename);
end

if isa(varargin{1},'wellObj') || isa(varargin{1},'MNW1Obj')
    well=varargin{1};
    varargin(1)=[];
else
    error('mfLab:plotWells:wrongObjectType',...
        'You must specify wells of class wellObj or MNW1Obj.');
end

if mustPlot
    for iw=length(well):-1:1

        if isempty(varargin)
            well(iw).whdl(1) = plot3(ax,well(iw).x,well(iw).y,well(iw).z,'color',defaultcolor,'linewidth',defaultlinewidth);
            well(iw).whdl(2) = plot3(ax,well(iw).x([1 1]),well(iw).y([1 1]),[well(iw).ztop well(iw).z(1)],'color',grey,'linewidth',1,varargin{:});
            well(iw).whdl(3) = plot(ax,well(iw).x,well(iw).y,'ko','markerfacecolor',defaultcolor);
            well(iw).UserData.color           = defaultcolor;
            well(iw).UserData.linewidth       = defaultlinewidth;
            well(iw).UserData.markerFaceColor = defaultcolor;
        else
            well(iw).whdl(1) = plot3(ax,well(iw).x([1 1]),well(iw).y([1 1]),well(iw).z([1 end]),varargin{:});
            well(iw).whdl(2) = plot3(ax,well(iw).x([1 1]),well(iw).y([1 1]),[well(iw).ztop well(iw).z(1)],'color',grey,'linewidth',1,varargin{:});
            i = strmatchi('color',varargin);
            if i,
                well(iw).UserData.color = varargin{i+1};
                well(iw).UserData.markerFaceColor = well(iw).UserData.color;
            else
                well(iw).UserData.color = defaultcolor;
                well(iw).UserData.markerFaceColor = well(iw).UserData.color;
            end
            well(iw).whdl(3) = plot(ax,well(iw).x(1),well(iw).y(1),'ko','markerfacecolor',well(iw).UserData.markerFaceColor);

            
            i = strmatchi('linewidth',varargin);
            if i,
                well(iw).UserData.linewidth = varargin{i+1};
            else
                well(iw).UserData.linewidth = defaultlinewidth;
            end
        end
    end
else
    for iw=1:length(well)
        if isfield('color',well(iw).UserData)
            color=well(iw).UserData.color;
        else
            color = grey;
        end
        
        if isfield('linewidth',well(iw).UserData)
            linewidth = well(iw).UserData.linewidth;
        else
            linewidth = 1;
        end
        
        if isfield('markerFaceColor',well(iw).UserData)
            markerFaceColor = well(iw).UserData.markerFaceColor;
        else
            markerFaceColor = grey;
        end
        
        set(well(iw).whdl(1),'visible','on');
        set(well(iw).whdl(2),'visible','on');
        if isnan(well(iw).Q(it))
           set(well(iw).whdl,'visible','off');
           
        elseif well(iw).Q(it)==0
           set(well(iw).whdl(1),'color',color);
           set(well(iw).whdl(1),'linewidth',linewidth);
           set(well(iw).whdl(2),'color',color);
           set(well(iw).whdl(2),'linewidth',linewidth);
           set(well(iw).whdl(3),'markerfacecolor',markerFaceColor);
        end
        
    end
end

end
