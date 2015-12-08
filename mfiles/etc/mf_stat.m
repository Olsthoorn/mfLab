function [m,M,S,iMin,iMax]=mf_stat(X,n)
%MF_STAT statistics of input X
%
% USAGE:
%     [min Max Std]=stat(X,n)
%
% X an array
% n is direction of dimension
%
% TO 100526 111220
%
% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if isnumeric(X)
    switch nargin
        case 1
            m=min(X(:));    
            M=max(X(:));
            S=std(X(:));
        case 2
            m=min(X,[],n);
            M=max(X,[],n);
            S=std(X,[],n);
        otherwise
            error('%s must have at least one argument',mfilename);
    end
    return;
end

if isstruct(X) && isfield(X,'values')
    m=Inf; M=-Inf; iMin=1; iMax=1;
    for it=1:numel(X)
        MIN = min(X(it).values(:));
        MAX = max(X(it).values(:));
        S=[];
        
        if MIN<m,
            iMin = it;
            m    = MIN;
        end
        
        if MAX>M
            iMax = it;
            M    = MAX;
        end
    end
end