%% For the new khettara heights
% Create a circle around the ending points of the khettaras (and also the
% outcrop). Interpolate the values of the ground surface. This is used to
% determine the ending height of the khettaras. The DEM map is showning to
% much uncertainty of the ground surface, and therefor the khettara height.

if 0
    %% This part has been done, and is saved to excel!!
    figure; hold on;
    [c,h] = gr.contourf(gr.Z(:,:,1),min(gr.Z(:)):1:850);
    title('Choose your circle'); axis equal;
    xlabel('x UTM [m]'); ylabel('y UTM [m]');
    hb = colorbar; set(get(hb,'title'),'string','Phi [m]');
    khettaras.plot
    modelArea.plot
    outcrops.plot

    [xx,yy] = ginput;
    [x,y] = utm2wgs(xx,yy,'30N');
    Surface = interp2(gr.Xm,gr.Ym,gr.Z(:,:,1),xx,yy);
    ExcelInfo = [x(:) y(:) xx(:) yy(:) Surface(:)];
end

figure; hold on;
Circle = area2Obj(gr,[XLSData, 'JorfData.xls'],'CircleKhet');
z = zeros(size(gr.Xm));
for i = 1:length(Circle.cellAreaVals)
    xU = Circle.Idx(i);
    z(xU) = Circle.cellAreaVals(i,5);
end
gr.contourf(z,min(gr.Z(:)):1:850);
hb = colorbar; set(get(hb,'title'),'string','Surface [m]');
khettaras.plot


for iK = 1:length(khettaras)
    if any(khettaras(iK).P(end).idx == Circle.Idx)
        display(khettaras(iK).name);
        z(iK) = Circle.cellAreaVals(find(khettaras(iK).P(end).idx == Circle.Idx),5);
        display(z(iK));
    else
        display(khettaras(iK).name);
        z(iK) = 0;
        display(z(iK));
    end
end