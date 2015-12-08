function [o,L]=plotLine(o,well)
% wellSeries = wellSeries.plotWells(well) --- plots wellSeries each in its own
% color as 'o' on XY plane
%
% L is lenght of each well series
%
% TO 120512

h=NaN(size(o));
L=NaN(size(o));

for is=1:length(o)
    I = find([well.parent]==o(is).id);
    if ~isempty(I)
        h(is)=plot([well(I).x],[well(I).y],[mf_color(is) '-']);
        L(is)=Length([well(I).x],[well(I).y]);
    end
end

end

function L = Length(X,Y)
    dx=diff(X(:));
    dy=diff(Y(:));
    L = sum(sqrt(dx.^2+dy.^2));
end