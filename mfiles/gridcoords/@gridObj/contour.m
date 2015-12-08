function [c,h] = contour(o,varargin)
%GRIDOBJ/SHOW: contours given array of proper size
%
% USAGE: h = gr.contourf(A,varargin)
%
% SEE ALSO:  gr.contour gr.contourf gr.contourXS, gr.contourYS gr.contourfXS gr.contourfY

[ax,varargin] = getNext(varargin,'axis',gca);
[A ,varargin] = getNext(varargin,'double',[]);

if ~isempty(A)
[c,h] = contour(ax,o.xc,o.yc,A,varargin{:});
end