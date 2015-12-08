function well = plotWells2D(well,varargin)
% well = plotWells2D(well,ax)     plots or updates the wells for SP sp
% well = plotWells2D(well,-iper)
%
% to refresh the wells, i.e to fully color the screen use -iPer (notice the
% minus sign). The minus sign is required to differentiate between a handle
% and the stress period number. This is only a small problem.
% When the well flow is not zero, the screen will be clored (white unless
% you change it by extra options) and when the well is off the screen is
% transparent with only its outline visible.
% You can provide any additional propertie/value combination int the
% argument list to change the properties of patch that represents the
% screen and the casing. (You can only access the two separately by
% changing it direcly in the well object. i.e. by a statement like
% set( well.whdl(1), property, value, property, value (for the casing)
% set( well.whdl(2), property, value, property, value (foor the screen)
%
% TO 120422

for iw=1:length(well)
    well(iw) = well(iw).plot2D(varargin{:});
end
