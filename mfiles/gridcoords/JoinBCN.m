function L=JoinBCN(L,JoinArray,code)
%JOINBCN joins a MODFLOW or MT3DMS/SEAWAT boundary condition list L of form [iPer iLay iRow iCol rest]
%
% Example:
%     L=JoinBCN(L,JoinArray,code)
%
%     used for refining MODFLOW/MT3D  models
%     Joining is according to the new layers implied by JoinArray
%
%     Nz is number of layers in old model
%     JoinArray has form [1:NOld; N1 N2 N2 N3 ...]
%     where N1 N2 are the net layer numbers (Ni<NOld)
%
%     As the new set of layers is a subset of teh old set, we only need to
%     replace the indices in the old array for the L, R or C.
%
% SEE ALSO: gridObj modelObj cutoutBCN cutoutXTA RefineBCN RefineGrid JoinBCN JoinLayers removeCBD
%
%  TO 100601 100610

switch lower(code(1))
    case 'x', LCol=4;
    case 'y', LCol=3;
    case 'z', LCol=2;
    otherwise
        error('Code argument in JoinBCN must be ''x'', ''y'' or ''z''.');
end

NOld=size(JoinArray,2);
for iOld=1:NOld  % length(NewLayers) is old Nz
    iNew=JoinArray(2,iOld);
    I=find(L(:,LCol)==iOld);
    if ~isempty(I)
        L(I,LCol)=iNew;
    end
end
