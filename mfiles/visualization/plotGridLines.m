function plotGridLines(varargin)
%PLOTGRIDLINES -- plot gid lines in an axes
%
% USAGE n = plotGridLines(ax,clr,alpha)

[ax   ,varargin] = getProp(varargin,'axis',gca);
[clr  ,~       ] = getNext(varargin,{'char,double'},grey);

xLim = get(ax,'xLim');
yLim = get(ax,'yLim');

if strcmpi('log',get(ax,'xScale'))
    xL = floor(log10(xLim(  1)));
    xR = ceil( log10(xLim(end)));
    xTick = (1:10)' * 10.^(xL:xR-1);
    xTick=unique(xTick(:));
    xTick=xTick(xTick>=xLim(1) & xTick<=xLim(end));
else
    xTick = get(ax,'xTick');
end

for i=1:numel(xTick)
    plot(ax,xTick([i i]),yLim,'color',clr); %,'edgeAlpha',alpha);
end

if strcmpi('log',get(ax,'yScale'))
    yB = floor(log10(yLim(  1)));
    yT = ceil( log10(yLim(end)));
    yTick = (1:10)' * 10.^(yB:yT-1);
    yTick=unique(yTick(:)); yTick=yTick(yTick>=yLim(1) & yTick<=yLim(end));
else
    yTick = get(ax,'yTick');
end

for i=1:numel(yTick)
    plot(ax,xLim,yTick([i i]),'color',clr); %,'edgeAlpha',alpha);
end
