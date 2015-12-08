function [c,h] = ge_contour(o,varargin)
%GRIDOBJ/SHOW: contours given array of proper size usis googleearth toolbox
%
% USAGE: h = gr.contourf(A,varargin)

[ax,varargin] = getNext(varargin,'axis',gca);
[A ,varargin] = getNext(varargin,'double',[]);

if ~isempty(A)
[c,h] = ge_contour(ax,o.ec,o.nc,A,varargin{:});
end