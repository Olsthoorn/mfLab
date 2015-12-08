function L = mf_find(A,which,dim)
%MF_FIND finds start or end of a zone in a multidimensional array along dimension dim
%
% Usage:
%    L = mf_find(A,which,dim)
%
%    A        a zone array with only true (or 1) and false (or 0)
%    L        is logical array of size A indicating first of last cells
%    which   {'first' | 'last} 
%    dim     { 1 | 2  | 3 }  dimension along which to search
%
% Example:
%     Let A be a zone array indicating some condition like WET>0 or IBOUND==1
%     so that A only contains ones (true) where the condition is met and
%     zeros (false) otherwise
%     Then to find the beginning or end of the zone in any dimension
%     LastZ  = mf_find(A,'last', 3);
%     FirstY = mf_find(A,'first',1);
%
% See also: find
%
% TO 120804

switch which
    case 'first'
        switch dim
            case 1
                L = diff(cat(dim,zeros(size(A(1,:,:))),A),1,dim);
            case 2
                L = diff(cat(dim,zeros(size(A(:,1,:))),A),1,dim);
            case 3
                L = diff(cat(dim,zeros(size(A(:,:,1))),A),1,dim);
            otherwise        
                error('%s: dim must be 1,2 or 3',mfilename);
        end
        L = (L==1);
    case 'last'
        A=flipdim(A,dim);
        switch dim
            case 1
                L = diff(cat(dim,zeros(size(A(1,:,:))),A),1,dim);
            case 2
                L = diff(cat(dim,zeros(size(A(:,1,:))),A),1,dim);
            case 3
                L = diff(cat(dim,zeros(size(A(:,:,1))),A),1,dim);
            otherwise        
                error('%s: dim must be 1,2 or 3',mfilename);
        end
        L=(L==1);
        L=flipdim(L,dim);        
    otherwise
        error('%s: which must be either first or last',mfilename);
end
