function list=cell2list(cellArray)
%CELL2LIST make list from a cell array of lists
%
% USAGE:
%    pntsrc=makelist(PNTSRC)
%
%    used in writeBCN to transfer the stresses as they have been collected
%    in an array of cells each with a list of modeflow-cell input. The cell
%    array is converted to a single sorted list that can be written to file
%    by the writefunction.
%
% TO 120410

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

    cellArray(cellfun(@isempty,cellArray)) = [];

    iHead = 5; % head column in most stress lists for modflow
    
    N=sum(cellfun('size',cellArray,1));
    W=max(cellfun('size',cellArray,2));
    
    list=NaN(N,W);
    
    k=1; I=find(~cellfun('isempty',cellArray));
    for i=I(:)'
        K=find(~isnan(cellArray{i}(:,min(W,iHead)))); % check head column
        if ~isempty(K)
            list(k:k+length(K)-1,:)=cellArray{i}(K,:);
            k=k+length(K);
        end
    end
    list=sortrows(list(~isnan(list(:,min(W,iHead))),:));
end
