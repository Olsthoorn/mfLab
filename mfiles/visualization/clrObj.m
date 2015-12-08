classdef clrObj
%CLROBJ class definition for color objects to plot on multiple axis with distinct colors (obsolete)
%
% Important to set colors and allowing to plot on multiple axes for
% sophisticated multi-purpose figures simultaneously
%
% See also: mf_setmulticolormap
%
% TO 110810
    properties
        ax    % axis handle
        range % data range for this axis
        crnge % contours levels to draw
        map   % colormap
        hdl   % handle of output of contour function
        fhdl  % handle to contour functiion used
        
        x     % plot axis
        y     % plot axis
    end
    methods
        function obj=clrObj(varargin)
            % obj=clrObj([axclr,]x,y,range,crnge,map)
            switch nargin
                case 0, return
                case 5 % backward compatibility
                    axclr     = 'none';
                case 6 % normal numer of arguments, do nothing
                    axclr    = varargin{1};
                    varargin = varargin(2:end);
                otherwise
                    error('%s: wrong number of input arguments (%d).',...
                        class(obj),nargin);
            end
            obj.ax    = axes('color',axclr); hold on;
            obj.x     = varargin{1};
            obj.y     = varargin{2};
            obj.range = varargin{3};
            obj.crnge = varargin{4};
            obj.map   = varargin{5};
        end
        function obj=contour(obj,funhdl,values)
            obj.fhdl=funhdl;
            
            [~,obj.hdl]=obj.fhdl(obj.ax,obj.x,obj.y,values,obj.crnge);
            if strcmp(func2str(obj.fhdl),'contourf')
                set(get(obj.hdl,'children'),'edgecolor','none');
            end
        end
        function obj=update(obj,values)
            set(obj.hdl,'zdata',values);
            if strcmp(func2str(obj.fhdl),'contourf')
                set(get(obj.hdl,'children'),'edgecolor','none');
            end
        end
        function obj=title(obj,ttl)
            title(obj.ax,ttl);
        end
    end
end