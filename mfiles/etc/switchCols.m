function A = switchCols(A,I,J)
% switch columns I and J in 2D array
% TO 130819

if numel(size(A))~=2
    error('%s: numel(size(A))  must be 2',mfilename,numel(size(A)));
end

if numel(I) ~=  numel(J)
    error('%s: Index vectors with number of columns to switch must be of equal length',mfilename);
end

if any(I)<1 || any(I)>size(A,2)
    error('%s: index1 must be >0 and < %d',mfilename,size(A,2));
end
if any(J)<1 || any(J)>size(A,2)
    error('%s: index2 must be >0 and < %d',mfilename,size(A,2));
end

B = A(:,I); A(:,J) = A(:,I); A(:,J) = B;
