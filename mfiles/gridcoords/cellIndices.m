function indices=cellIndices(I,dims,orderstr)
%CELLINDICES get individual axes indices given global one, chosen order
%
% USAGE:
%   indices=cellIndices(I,dims,orderstr)
%
%   Get the indices in an n-dimensional grid given by its dims=size(M)
%   orderstr is a string to change the order of the columns using L R C
%   for layer (initially col 3), row (initially col 1) and C (initially col 2)
%   if number of dimensions through dims <>3, then orderstr is ignored.
%
% EXAMPLE:
%    LRC = cellIndices(I,size(A),'LRC');
%    RLC = cellIndices(I,size(A),'RLC');
%    etc.
%   LRC is a list of the nodes [Layer Row Col] A is a 3D array
%   and I are the global indices in the array A
%
% SEE ALSO: also cellIndex xyzindex linegrid inpolygon IdxMatlab2Modflow
%
% TO 090317 091119


% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if length(dims)<length(orderstr), dims(end+1:length(orderstr))=1; end
    
if islogical(I), I=find(I); end

indices=zeros(length(I),length(dims));
for i=1:length(I)
    n=I(i)-1;
    for idim=1:length(dims)
        indices(i,idim)=rem(n,dims(idim))+1;
        n=(n-(indices(i,idim)-1))/dims(idim);
    end
end

if exist('orderstr','var')
    if length(orderstr)==3
        switch upper(orderstr)
            case 'LRC', order=[3 1 2];
            case 'LCR', order=[3 2 1];
            case 'RLC', order=[1 3 2];
            case 'RCL', order=[1 2 3];
            case 'CLR', order=[2 3 1];
            case 'CRL', order=[2 1 3];
            otherwise error(sprintf('Order string <<%s>> in cellIndices',orderstr));
        end
    elseif length(orderstr)==2
        switch upper(orderstr)
            case 'RC', order=[1 2];
            case 'CR', order=[2 1];
            otherwise error(sprintf('Order string <<%s>> in cellIndices',orderstr));
        end
    else
        error('orderstr must be a combination of the letters LRC or RC');
    end
    indices=indices(:, order);  % switch colmns
end
