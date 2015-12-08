function o = concYS(o,varargin)
%ANIMATEOBJ/CONCYS -- shows transient conc in cross section along y-axis
%
% USAGE:
%    animateObj = animateObj.concXS(gr,IComp,I,well,STCONC,varargin)
%
% gr       = gridObj
% IComp    = index/Nr of compoments/Species to plot
% I        = rows (average conc of the selected rows will be used)
% well     = wells or MNW (multi node wells, of class MNW1Obj)
% STCONC   = STCONC may be omitted
% varargin = parameters to pass to the plotting routines and options
% options that can be used as option,value pairs:
%    'figure' figname
%    'figPos' figpos [xLL yLL w h] in pixels
%    'STCONC' STCONC
%    'xlim' xlim
%    'zlim' zlim (vertical dimension to show)
%    'VK',VK to compute and visualize the hydraulic resistance of layers.
%
% if options are omitted --> defaults will be used:
%  figure = gcf
%  figPos = screenPos(0.75), that is 75% of the screen size
%  STCONC = [];
%  xLim   = gr.xc([1 end]);  same as gr.xGr([1 end])
%  yLim   = gr.zc([end  1]); same as gr.zGr([end 1])
%  I   = rows of wells to be shown, if omitted and if well objects exist
%  I   = 1 if model has only one row
%  IComp = all components 1:NCOMP
%
% TO 120908 121017 130403

o = o.conc2D('y',varargin{:});

