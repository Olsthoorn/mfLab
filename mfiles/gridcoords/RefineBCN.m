function BCNnew=RefineBCN(BCNold,SplitArray,code)
%REFINEBCN refines stresses lists WEL, DRN, RIV, CHD according to SplitArray
%
% Example:
%    BCNnew = RefineBCN(BCNold,SplitArray,code)

%    Refines WEL, DRN, RIV, CHD according to SplitArray and code,
%    where code is 'x', 'y' or 'z' and defined the direction to split
%    SplitArray has the following form [i1 i2 i3 i4; N1 N2 N3 N4; I1 I2 I3 I4]
%    The first line are the cell indices of the current model to be
%    subdivided (splitted up). Only the cells to be splitted need to be
%    mentioned.
%
%    The second line indicates into how many new cells this old cell will be
%    splitted. Any positive integer is acceptable.
%
%    The third line is optoinal.
%    If it is ommitted, then the conductances of the boundary list
%    or the flow of a well list will be equally divided over the new refined
%    grid cells.
%
%    Else it will indicate over which new cells the conductances and the flows
%    of the well list will be (equally) distributed.
%    0  = the center of the new subcells gets the complete concuctance or flow
%    <0 = The conductance or flow is equally divided over the first new cells
%       up to the number given by the value. If this number is larger than the
%       the number of subcells indicated in the second line it is the same as
%       division equally over all cells.
%    >0 = The conductances or flow is equally divided over the last new cells
%       as indicated by the number given on the third line. Use -1 to
%       attribute all weight on the last cell, -4 to divie the total
%       conductance or flow over the last 4 of thenew cells. If abs(value)
%       is larger than the total number of new cells, then the conductance
%       or flow will be equally distributed over all new cells. The result
%       will then be the same as a large positive number > NNew or omitting
%       the third line of Splitarray altogether.
%
% These function are superseded by methods of the gridObj
%
% SEE ALSO: gridObj modelObj cutoutBCN cutoutXTA RefineBCN RefineGrid JoinBCN JoinLayers removeCBD
%
%
% TO 100604 120514

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

code=lower(code(1));

switch lower(code(1))
    case 'x', Lcol=4;
    case 'y', Lcol=3;
    case 'z', Lcol=2;
    otherwise
        error('code must be ''x'',''y'' or ''z'' to indicate splitting direction');
end

Nold = max(max(SplitArray(1,:)),max(BCNold(:,Lcol)));

oldCells2Split = SplitArray(1,:);
nrOfSubCells   = SplitArray(2,:);
subCellPos     = SplitArray(3,:);

%% Grab the BCN in the column (no matter which y or z
%  Get list of existing indices for every column

Idxold(Nold)=struct('Connected',[],'NSplit',[],'weight',[],'iNew',[]); % Struct to store BCNold info

N   =0; % Counts required number of BCNnew for alloccation of mem
iNew=0; % Index of new BCN to replace values in Lcol

for iCell=1:Nold % number of cells along the axis to be refined
    
    Idxold(iCell).Connected=find(BCNold(:,Lcol)==iCell);  % BCNolf connected to this Index
    
    j = find(oldCells2Split==iCell);
    
    if ~isempty(j)
        Idxold(iCell).NSplit = nrOfSubCells(j);
        
        % New cell numbers linked to old cell i:
        Idxold(iCell).iNew   = iNew + (1:Idxold(iCell).NSplit);
        
        Idxold(iCell).weight=zeros(1,Idxold(iCell).NSplit);
        if subCellPos(j)<0     % distribute BCNold over the first abs(subCellPos(j))
            Idxold(iCell).weight(1:abs(subCellPos(j)))=1/abs(subCellPos(j));
        elseif subCellPos(j)>0 % distributed BCNold over the last abs(subCellPos(j)
            Idxold(iCell).weight(end-subCellPos(j)+1:end)=1/subCellPos(j);
        else % put total BCNold in center of BCNnew
            imid = max(1,floor(Idxold(iCell).NSplit/2));
            Idxold(iCell).weight(imid) = 1;
        end
        
         iNew=iNew + Idxold(iCell).NSplit;
    else % cell not split
        Idxold(iCell).NSplit = 1;
        Idxold(iCell).iNew   = iNew+1;
        Idxold(iCell).weight = 1;
        
         iNew=iNew+1;
    end
    N=N+sum(Idxold(iCell).weight>0)*length(Idxold(iCell).Connected);
end

BCNnew = NaN(N,size(BCNold,2));

iCOND = 6; % The column holding the conductance

k=0;
for iCell=1:length(Idxold)
    if ~isempty(Idxold(iCell).Connected)
        for j=find(Idxold(iCell).weight>0)
            BCNnew(k+(1:size(Idxold(iCell).Connected,1)),   :) = BCNold(Idxold(iCell).Connected,:);
            BCNnew(k+(1:size(Idxold(iCell).Connected,1)),Lcol) = Idxold(iCell).iNew(j);
            if ~code=='z'
                BCNnew(k+(1:size(Idxold(iCell).Connected,1)),iCOND) = BCNold(Idxold(iCell).Connected,iCONC) ...
                                                                     *BCNold(Idxold(iCell).weight(j));
            end
            k =    k+   size(Idxold(iCell).Connected,1);
        end
    end
end
