function o=plotWells(o,well)
% wellSeries = wellSeries.plotWells(well) --- plots wellSeries each in its own
% color as 'o' on XY plane
%
% TO 120512

h=NaN(size(o));
for is=1:length(o)
    I = find([well.parent]==o(is).id);
        
    fprintf('wellSeries %2d, Nr wells belonging %d\n',o(is).id,sum([well.parent]==o(is).id));
    
    if ~isempty(I)
        h(is)=plot([well(I).x],[well(I).y],[mf_color(is) 'o']);
    end
end