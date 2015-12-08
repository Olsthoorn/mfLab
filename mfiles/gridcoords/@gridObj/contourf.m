function [c,h] = contourf(o,varargin)
%GRIDOBJ/SHOW: contours given array of proper size
%
% USAGE: h = gr.contourf(A,varargin)
%
% SEE ALSO: gr.contourXS, gr.contourYS gr.contourf gr.contourfXS gr.contourfY

[ax,varargin] = getNext(varargin,'axis',gca);
[A ,varargin] = getNext(varargin,'double',[]);

if ~isempty(A)
[c,h] = contourf(ax,o.xc,o.yc,A,varargin{:});
end