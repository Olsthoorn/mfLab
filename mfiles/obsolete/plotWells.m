function h=plotWells(~,varargin)
%% h = plotWells(ax,wells)  -- plots wellscreens in 3D
%
% TO 120531

if ~exist('varargin','var') || isempty(varargin)
    error('gridObj:plotWells:noData',...
        'gridObj/%s: argument wellObj needed',mfilename);
end

if ~isaxis(varargin{1})
    ax=gca;
else
    ax = varargin{1};
end
varargin(1)=[];

if isempty(varargin)
    error('grid:plotWells:noWellData',...
        'gridObj/%s: argument wellObj needed',mfilename);
end
well=varargin{1}; varargin(1)=[];

if ~(isa(well,'wellObj') || isa(well,'MNW1Obj')
    error('mfLab:plotWells:wrongObjectType',...
        'You must specify wells of class wellObj.');
end

set(ax,'nextplot','add');

h = NaN(size(well));

for iw=1:length(well)

    if isempty(varargin)
        h(iw)=plot(well(iw).x([1 1]),well(iw).y([1 1]),well(iw).z,'b','linewidth',3);
    else
        h(iw)=plot(well(iw).x([1 1]),well(iw).y([1 1]),well(iw).z,varargin{:});
    end
end

end
