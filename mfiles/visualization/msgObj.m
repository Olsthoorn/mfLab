classdef msgObj
    %% Places a message box on the figure and allows dynamic updating of the
    %% message in a loop.
    properties
        type='msgObj';
        hdl;
        x; y;
        string;
    end
    methods
        function o=msgObj(varargin)
            if nargin==0
                return;
            end
            if ishandle(varargin{1}),
                axes(varargin{1});
                varargin=varargin(2:end);
            end
            if length(varargin)<3
                error('mfLab:msgObj:notEnoughInputArguments',...
                    'msgObj: not enough input arguments: use {ax,x,y,txt} or {x,y,txt}\n');
            end
            o.x      = varargin{1};
            o.y      = varargin{2};
            o.string = varargin{3};
            varargin=varargin(4:end);
            o.hdl=text(o.x,o.y,o.string,varargin{:});
            set(o.hdl,'fontname','courier');
        end
        function hdl=update(o,str)
            o.string=str;
            set(o.hdl,'string',str);
            hdl=o.hdl;
        end
    end
end


%% Text properties
% get(hdl)
% 	Annotation = [ (1 by 1) hg.Annotation array]
% 	BackgroundColor = none
% 	Color = [0 0 0]
% 	DisplayName = 
% 	EdgeColor = none
% 	Editing = off
% 	Extent = [0.995402 0.973761 0.0390805 0.0408163]
% 	FontAngle = normal
% 	FontName = Helvetica
% 	FontSize = [10]
% 	FontUnits = points
% 	FontWeight = normal
% 	HorizontalAlignment = left
% 	LineStyle = -
% 	LineWidth = [0.5]
% 	Margin = [2]
% 	Position = [1 1 0]
% 	Rotation = [0]
% 	String = jan
% 	Units = data
% 	Interpreter = tex
% 	VerticalAlignment = middle
% 
% 	BeingDeleted = off
% 	ButtonDownFcn = 
% 	Children = []
% 	Clipping = off
% 	CreateFcn = 
% 	DeleteFcn = 
% 	BusyAction = queue
% 	HandleVisibility = on
% 	HitTest = on
% 	Interruptible = on
% 	Parent = [178.051]
% 	Selected = off
% 	SelectionHighlight = on
% 	Tag = 
% 	Type = text
% 	UIContextMenu = []
% 	UserData = []
% 	Visible = on