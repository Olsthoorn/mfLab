function stress = mvStressDown(o,stressName,stress,minThickness)
% stress = stressDown(stressName,stress,minThickness) - moves
%   stress z-cell down if layer is too thin
%   This proved necessary in NHI submodels where top layers are
%   convertible and often very thing, preventing convergence. If the stress
%   is put in a lower (thicker) layer by increasing z index of the stress
%   by 1, the convergence issue was solved.
%
% NHI is the Dutch National Hydrologic Instrument, a nation-wide
% groundwater-surfacewater model,whose groundwater model data are published
% on www.NHI.nu for the grids at a resolution of 250x250 m.
% see mfLab/mfiles/NHI and mfLab/examples/NHI
%
% TO 120816

if nargin<4, minThickness = 1.0; end

    fprintf('Increase layer number of stress %s if layer thickness < %g m\n',stressName,minThickness);

    fmt ='%5d %3s cells lowered to underlying layer.\n';
    % Get global indices of stress cells
    Idx = cellIndex(stress(:,[4 3 2]),o.size);
    if ~isnan(Idx)
        tooThin = o.DZlay(Idx)<minThickness;
        stress(tooThin,2)=stress(tooThin,2)+1;
        fprintf(fmt,sum(tooThin),stressName); 
    end    
