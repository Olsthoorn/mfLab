function [isLay,LAYCBD]=isLayer(Nz,LAYCBD)
% isLay = is a vector with 1 if a model layer and 0 if confining bed
% TO 120414

%% Clean up LAYCBD
LAYCBD(end:Nz)=0;  % make LAYCBD sufficiently long
LAYCBD = logical(LAYCBD(:)); % make sure only 0 and 1 are used
laycbd = [ones(size(LAYCBD)),LAYCBD]'; % complete with model layers as first column
laycbd = reshape(cumsum(laycbd(:)),[2,numel(LAYCBD)]); % count z-layers
LAYCBD(min(laycbd,[],1)>Nz)=[]; % cut where where z-layer>Nz

if LAYCBD(end)
    % last LAYCBD must be false, if not replace by layer
    LAYCBD= [LAYCBD(1:end-1); false ; false];
end
% cleaned

isLay=false(Nz,1); % true of z-layer is a model layer

iz=1;  % index of z-layer (not model layer)
for iLay=1:length(LAYCBD)
    if LAYCBD(iLay)
        isLay(iz)=true;  iz=iz+1; if iz>Nz, break; end
        isLay(iz)=false; iz=iz+1;
    else
        isLay(iz)=true; iz=iz+1;
    end
    if iz>Nz, break; end
end
