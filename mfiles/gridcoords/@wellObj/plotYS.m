function well = plotYS(well,varargin)
% well = well.plotYS(ax,ix,varargin)  plot wells on Y-section.
% well = well.plotYS(it)  refresh well status for this stress period.
% 
% TO 120510 120407

if isempty(varargin)
   well = well.plot2D(gca,'y');
   return;
else
    well = well.plot2D(varargin{1},'y',varargin{2:end});
end
