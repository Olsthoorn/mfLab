function L = sameSize(A,B)
%SAMESIZE check if size of A and B are the same
%
%USAGE: L = sameSize(A,B)
%   L is true of size of A and B are the same and A and B are of same class
%
% TO 130514

L = false;

% if ~strcmpi(class(A),class(B))
%     error('%s: A and B must be of same type',mfilename);
% end

if numel(size(A)) == numel(size(B))
    L =  all(size(A) == size(B));
end