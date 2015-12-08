function plotgrid(varargin)
%PLOTGRID plots the grid lines in color clr given the coordinates xGr yGr
%    and possibly well locations in blue if well is a struct
%    whose elements have fields x and y
%
% Example:
%    plotgrid([ax,]xGr,yGr,clr,lw,well,figname,figcoords)
%    plotgrid(xGr,yGr, [clr [,well [,finame,figcoords]])
%    plotgrid(xGr,yGr, clr,lw)  % with clr clr-line type combination and lw linewidth
%
%   if figname is given with figcoords, then fig is plotted too
%   if figcoords are missing xGr and yGr will be used
%
% ToDo: remove in favor of method plotGrid of gridObj
%
% See also: gridObj
%
% TO 091201 100115



try
    get(varargin{1},'xgrid','linear');
    set(axes,varargin{1});
    varargin(1) =[];
catch %#ok
    ax=gca;
end

xGr=varargin{1};
yGr=varargin{2};
varargin(1:2)=[];

if ischar(varargin{1})
    clr = varargin{1};
elseif numel(varargin)==3
    clr = varargin{3};
else
    clr='c';
end

i = strmatchi('linewidth',varargin);
if i
    lw =varargin{i+1};
else
    lw =1;
end

i = strmatchi('well',varargin);
if i
    well = varargin{i+1};
end

i = strmatchi('figco',varargin);
if i
    figcoords = varargin{i+1};
    xgridpic= figcoords([1 3]);
    ygridpic= figcoords([2 4]);
else
    xgridpic= xGr([1 end]);
    ygridpic= yGr([1 end]);
    ygridpic=sort(ygridpic);
end

i = strmatchi('figname',varargin);
if i
    figname = varargin{i+1};
    gridpic = imread(figname);
    image(xgridpic,ygridpic,gridpic);
end


if min(size(xGr))==1, xGr=xGr(:)'; end
if min(size(yGr))==1, yGr=yGr(:) ; end

for i=1:size(xGr,2)
    if min(size(xGr))==1
        plot(ax,xGr([i,i]),yGr([1,end]),clr,'linewidth',lw);
    else
        plot(ax,xGr(:,i),yGr(:,i),clr,'linewidth',lw);
    end
end

for i=1:size(yGr,1)
    if min(size(yGr))==1
        plot(ax,xGr([1, end]),yGr([i,i]),clr,'linewidth',lw);
    else
        plot(ax,xGr(i,:),yGr(i,:),clr,'linewidth',lw);
    end
end

if exist('well','var')
    for iW=1:length(well), plot(ax,well(iW).x,well(iW).y,'or','markerfacecolor','r'); end
end
