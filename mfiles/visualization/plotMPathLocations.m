function plotMPathLocations(mpathLOCATIONS,varargin)
%PLOTMPATHLOCATIONS plot  MPath locations made by gr.startLoc in 3D
%
% Example:
%    plotMPathLocations(mpathLocations);
%
% TO 121202

    xyz  = gr.relloc2model(mpathLOCATIONS);
    
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
    
    if ~isempty(varargin)
        I = strmatchi('linestyle',varargin);
        if I
            i=I(1);
            lstyles = varargin{i+1};
            varargin(i:i+1)=[];
        else 
            lstyles={'-' '--' ';' '-.'};
        end
    else
            lstyles={'-' '--' ';' '-.'};
    end

    if isempty(varargin), varargin = {'visible' ,'on'}; end
    
    for i=1:numel(xyz)
        plot3(xyz{i}(:,1),xyz{i}(:,2),xyz{i}(:,3),[pick(i,clrs), pick(i,lstyles)],varargin{:});
    end
