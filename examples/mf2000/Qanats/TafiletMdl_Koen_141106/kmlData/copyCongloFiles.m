D = dir('conglo*.kml');

layerNames = {'limons','graviers','calcaire','conglomerate'};
for i=1:numel(D)
    for j=1:numel(layerNames)-1
        restName = D(i).name(13:end);
        copyfile(D(i).name,[layerNames{j} restName]);
    end
end