function M=xsec(M,varargin)
%XSEC: turns 3D array so that 3rd and 1st diemsnions are switched
% easy viewing and contouring of cross sections
%
% Example:
%   M=xsec(M [,dim])
%
%   if dim is '1' or 'r' permute M to [L,C,R]
%   Contouring of a column:
%   if dim is '2' or 'c' permute M to [L,R,C]
%
% TO 090101 110420

% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if       length(varargin)<1, dim=1;
elseif   ischar(varargin{1}),
    dim=lower(varargin{1});
else
    dim=varargin{1};
end

switch dim
    case {1 'r'},     M=permute(M,[3 2 1]);
    case {2 'c'},     M=permute(M,[3,1,2]);
    otherwise
        error('illegal dimension: use 1, ''row'' or 2 or ''col'' for dim');
end
