function Anew=JoinLayers(Aold,JoinArray,DZ,code)
%JOINLAYERS joins layer array OldLayer according to JoinArray
%
% Example:
%   NewArray=JoinLayers(OldLayer,JoinArray,DZold,code)
%
%   JoinArray has form [OldLayers; NewLayers] Such that any old layers
%   will be joined to one new layer.
%   DZold is the array with the DZ of the original model
%   code is 'k', 'c' or 'z'
%
%   k joins according to transmissivities
%   c joins according to resistances
%     joins by attributing to the new joined layer (for instance IRCH or
%     RECH that have one value per Ix,Iy location.
%   z yields the new cell elevations, then OldLayer=Zold and NewArray=Znew
%   Z is the array with top and bottom elevations of cells
%
% These function are superseded by methods of the gridObj
%
% SEE ALSO: gridObj modelObj cutoutBCN cutoutXTA RefineBCN RefineGrid JoinBCN JoinLayers removeCBD
%
% TO 100601

if ~all(JoinArray(1,:)==sort(JoinArray(1,:)))
    error('JoinArray(1,:) must be 1: number of old layers');
end

NewLayers=unique(JoinArray(2,:));

if ~all((1:length(NewLayers)) == NewLayers)
    error('JoinArray(2,:) must contain all new layers 1:NNew<NNold');
end

% size(NewArray) = [NyOld,NxOld,NzNew]
if code=='i' % e.g. for IRCH Ny*Nx*Nper with layer indices to put recharge
    Anew=Aold;
else
    Anew=NaN([size(Aold(:,:,1)) length(NewLayers)]);
end

for iLay=NewLayers
    IZ=find(JoinArray(2,:)==iLay);  % all old layers to be joined into NewLayer(iLay)
    switch code
    case 'k' % join by adding k*dz
        Anew(:,:,iLay)=sum(Aold( :,:,IZ).*DZ(:,:,IZ),3)./sum(DZ(:,:,IZ),3);
    case 'c' % join by harmonic mean (adding inverse)  dz/k
        Anew(:,:,iLay)=sum(DZ(:,:,IZ),3)./sum(DZ(:,:,IZ)./Aold( :,:,IZ),3);
    case 'z' % join by adding the thicknesses of the layer
        if size(Aold,3)==size(DZ,3)+1  % cell elevations
            if iLay==1,
                Anew(:,:,1)=Aold(:,:,1);
            end
            Anew(:,:,iLay+1)=Aold(:,:,IZ(end)+1);
        else
            Anew(:,:,iLay)=sum(Aold(:,:,IZ),3);
        end
    case 'i' % select new layer
        Anew(Aold>=IZ(1) & Aold<=IZ(end))=iLay; % attribute to new layer index
    otherwise
        error('code must be ''k'' or ''c'' or ''i'' or ''z''');
    end
end
