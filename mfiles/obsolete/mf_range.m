function [mM]=mf_range(arange,fldnm)
%MF_RANGE maximum and minimum value in range
%   range may be double, cell-array or struct array
%
%EXAMPLE
%   [m,M]=mf_range(srange);
%   [m,M]=mf_range(srange,fldnm);
%
% TO 110423

if isnumeric(arange),
    m=min(arange(:)); M=max(arange(:));
elseif iscell(arange),
    m=Inf; M=-inf;
    for i=1:numel(arange)
        m=min(m,arange(:));
        M=max(M,arange(:));
    end
elseif isstruct(arange)
    m=Inf; M=-Inf;
    if nargin<2, fldnm='values'; end
    for i=1:numel(arange)
        m=min(m,arange(i).(fldnm));
        M=max(M,arange(i).(fldnm));
        while 1, m=min(m); if numel(m)==1, break; end, end
        while 1, M=max(M); if numel(M)==1, break; end, end
    end
end

mM=[m M];
