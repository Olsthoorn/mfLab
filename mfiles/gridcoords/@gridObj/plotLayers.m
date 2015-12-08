function h=plotLayers(o,varargin)
% h = o.plotLayers(ax,layers,C,varargin) -- Plots layers of grid
%
% TO 120529

%% Get axis
if nargin==1,
    ax=gca;
elseif isaxis(varargin{1})
    ax=varargin{1};
    varargin(1)=[];
else
    ax=gca;
end

%% Getlayers
if isempty(varargin)
    layers=1:o.Nlay;
elseif isnumeric(varargin{1}) || islogical(varargin{1})
    layers = varargin{1};
    varargin(1) = [];
    if logical(layers), layers=find(layers); end
else
    error('You must specify layers in second argument');
end

%% Get data to plot
if isempty(varargin) || ~(isnumeric(varargin{1}) || isstruct(varargin{1}))
    error('%s: You must provide data (a 3D data array or a struct with field values)',mfilename);
end

C=varargin{1};
varargin(1) =  [];

crange = ContourRange(C,50); clim=[min(crange) max(crange)];
set(ax,'clim',clim);

%% plot the layers as 3D planes
h=NaN(size(layers));

set(ax,'nextplot','add');

for iLay = 1:length(layers)
    if isempty(varargin)
        if isstruct(C)
            h(iLay) = surf(ax,o.XGr,o.YGr,o.ZBgrlay(:,:,iLay),C(end).values(:,:,iLay));
        else
            h(iLay) = surf(ax,o.XGr,o.YGr,o.ZBgrlay(:,:,iLay),C(:,:,iLay));
        end
    else
        if isstruct(C)
            h(iLay) = surf(ax,o.XGr,o.YGr,o.ZBgrlay(:,:,iLay),C(end).values(:,:,iLay),varargin{:});
        else
            h(iLay) = surf(ax,o.XGr,o.YGr,o.ZBgrlay(:,:,iLay),C(:,:,iLay),varargin{:});
        end
    end
end
drawnow;

end
