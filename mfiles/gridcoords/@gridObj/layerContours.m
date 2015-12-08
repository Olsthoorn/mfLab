function layerContours(o,varName,var,ttl,dimension)
% gr.layerContours(varName,var,ttl,dimension) --- contour variable for all layers
% 
% TO 120816
if nargin<5, dimension=''; end
if nargin<4, ttl=[]; else  ttl = [ttl ': ']; end

range= ContourRange(var(~isnan(var)),50);

for iLay=1:size(var,3)
    figure; hold on; xlabel('x [m]'); ylabel('y [m]');

    title(sprintf('%s%s layer %d, dvalue=%4g %s',ttl,varName,iLay,min(diff(range)),dimension));
    
    o.contourf(gca,var(:,:,iLay),range);
    colorbar;
end
