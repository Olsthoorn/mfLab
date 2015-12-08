function BCN = cutoutBCN(BCN,Ix,Iy,LAYCBDold,LAYCBDnew)
%CUTOUTBCN cuts out boundary WEL, GHB,CHD,DRN,RIV, ... to match coordinate indices Ix and Iy
%
% Example:
%    BCN = cutoutBCN(BCN,Ix,Iy,LAYCBDold,LAYCBDnew)
%
%
% works for modelObj
%
% These function are superseded by methods of the gridObj
%
% SEE ALSO: gridObj modelObj cutoutBCN cutoutXTA RefineBCN RefineGrid JoinBCN JoinLayers removeCBD
%
% TO 120507

BCN=BCN( BCN(:,3)>=min(Iy) & BCN(:,3)<=max(Iy),:);
BCN=BCN( BCN(:,4)>=min(Ix) & BCN(:,4)<=max(Ix),:);

BCN(:,3)=BCN(:,3)-Iy(1)+1;
BCN(:,4)=BCN(:,4)-Ix(1)+1;

if nargin>4
    Iout = convertLAYCBD(LAYCBDold,LAYCBDnew);

    BCN(:,2) = Iout(BCN(:,2));
end
