function Z = mf_interp2(X,Y,Z,Ix,Iy,varargin)
%MF_INTERP2 -- switch to interp1 if X or Y is one col or one row resp.
%
% USAGE: see interp2
%
% Catch when indeed the interp2 is an interp1 as happens with
% models that consist of one row or one columns
% TO 131123

if size(Y,1)==1
    Z = interp1(X,Z,Ix,varargin{:});
elseif size(X,2)==1
    Z = interp1(Y,Z,Iy,varargin{:});
else
    Z = interp2(X,Y,Z,Ix,Iy,varargin{:});
end