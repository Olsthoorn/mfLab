function o = headYS(o,varargin)
%ANIMATEOBJ/HEADYS -- shows transient heads in cross section along y-axis
%
% USAGE:
%    animateObj = animateObj.headXS(gr,I,well,STRTHD,varargin)
%
% gr       = gridObj
% IComp    = index/Nr of compoments/Species to plot
% I        = column number or 1 if omitted
% well     = wells or MNW (multi node wells, of class MNW1Obj)
% STRTHD   = STRTHD may be omitted (causes drawdown to be plotted if present)
% varargin = parameters to pass to the plotting routines and options
% options that can be used as option,value pairs:
%    'figure' figname
%    'figPos' figpos [xLL yLL w h] in pixels
%    'STRTHD' STRTHD
%    'xlim' xlim
%    'zlim' zlim (vertical dimension to show)
%    'VK',VK to compute and visualize the hydraulic resistance of layers.
%
% if options are omitted --> defaults will be used:
%  figure = gcf
%  figPos = screenPos(0.75), that is 75% of the screen size
%  STRTHD = [];
%  xLim   = gr.xc([1 end]);  same as gr.xGr([1 end])
%  yLim   = gr.zc([end  1]); same as gr.zGr([end 1])
%  I   = columns of wells to be shown, if omitted and if well objects exist
%  I   = 1 if model has only one row
%
% TO 120908 121017 130403

o = o.head2D('y',varargin{:});

