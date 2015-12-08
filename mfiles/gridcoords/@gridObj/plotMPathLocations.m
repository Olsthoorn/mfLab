function plotMPathLocations(gr,mpathLocations,varargin)
% gridObj.plotMPathLocations(mpathLocations,varargin) --- 3D plot of MPaht locations made by gr.startLoc
% varargin -- e.g. 'color'.'rgbk',linestyle',{'--' '-' '_.'},lineWidth,[1 2 2 1]
% TO 121202

    xyz  = gr.relloc2model(mpathLocations);
    
    %% Color
    if ~isempty(varargin)
        I = strmatchi('color',varargin);
        if I
            i=I(1);
            clrs = varargin{i+1};
            varargin(i:i+1)=[];
        else
            clrs = 'brgk';
        end
    else
        clrs = 'brgk';
    end
    
    %% Linestype
    if ~isempty(varargin)
        I = strmatchi('marker',varargin);
        if I
            i=I(1);
            marker = varargin{i+1};
            varargin(i:i+1)=[];
        else 
            marker={'.' 'o' 'x' '+' 's' 'p' '^' 'v'};
        end
    else
        marker={'.' 'o' 'x' '+' 's' 'p' '^' 'v'};
    end
    
    %% Linewidth
    if ~isempty(varargin)
        I = strmatchi('linewidth',varargin);
        if I
            i=I(1);
            markersize = varargin{i+1};
            varargin(i:i+1)=[];
        else 
            markersize= 5;
        end
    else
        markersize=5;
    end

    if isempty(varargin), varargin = {'visible' ,'on'}; end
    
    for i=1:numel(xyz)
        plot3(xyz{i}(:,1),xyz{i}(:,2),xyz{i}(:,3),[pick(i,clrs), pick(i,marker)],...
            'markersize',pick(i,markersize),varargin{:});
    end
