function well=plotXS(well,varargin)
% well = well.plotXS(ax,ix,varargin)  plot wells on X-section.
% well = well.plotXS(it)  refresh well status for this stress period.
% 
% TO 120510 120407

well = well.plot2D('x',varargin{:});
