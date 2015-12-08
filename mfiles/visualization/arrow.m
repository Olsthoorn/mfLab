function hdl = arrow(varargin)
%ARROW draw/plot an arrow
%
% Example:
%    hdl= arrow(ax,x,y,npx,angle,varargin)
%
% TO 120505

[ax   ,varargin] = getNext(varargin,'axis',gca);
[x    ,varargin] = getNext(varargin,'double',0);
[y    ,varargin] = getNext(varargin,'double',0);
[npx  ,varargin] = getNext(varargin,'double',8);
[angle,varargin] = getNext(varargin,'double',8);

%% Arrow size and shape
if isempty(varargin), varargin = {'color','b','linewidth',2}; end

%% Arrow in pixels, without rotation
dxpx= npx * [ 0  5  5  8  5  5  0];
dypx= npx * [-1 -1 -3  0  3  1  1];

%% Arrow in pixels, rotated (angle is in as viewed by user on screen anti clockwise)
[dxpx,dypx] = rotate(dxpx,dypx,dxpx(end),dypx(end),angle);

%% Put the arrow on the screen
axUnits = get(ax,'Units');  % remember user units
set(ax,'Units','pixels');   % change to pixel units

%% limits in user coordinates
xlim=get(ax,'xlim');
ylim=get(ax,'ylim');

%% arrow in pixels on screen
figPos = get(gcf,'position'); % get figure position (always in pixel units)
axPos  = get(ax ,'position');

fx = (xlim(2)-xlim(1))/axPos(3);  % axPos(3) in pixels
fy = (ylim(2)-ylim(1))/axPos(4);  % axPos(4) in pixels

hdl = NaN(size(x));

[color,varargin] = getProp(varargin,'color','k');
[color,varargin] = getNext(varargin,'char',color);
[color,varargin] = getNext(varargin,'double',color);
for i=numel(x)
    hdl(i) = fill(x(i)+fx*dxpx, y(i)+fy*dypx,color,varargin{:});
end

% restore axis units
set(ax ,'Units',axUnits);