function well = plotYS(o,varargin)
% gr.plotYS(ax,ix,well,varargin)  plot wells on Y-section.
% gr.plotYS(it,well)  refresh well status for this stress period.
% 
% TO 120510 120407

error('%s: TO121119 THis method is under the wrong class, must be with class wellOb not gridObj',mfilename);


well=o.plot2D('y',varargin{:}); % must plot
