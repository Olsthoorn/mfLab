function X=mf_expand(IV,dims)
%MF_EXPAND expands matrix to full matrix
%
% USAGE:
%    X = mf_expand([I, V],dims)
%
%    Expands matrix to full matrix. Where IV is matrix given in compact form
%    i.e. as [I V] where both I and V are column vectors with I the global
%    coordinate indices and V the actual values. The rest of the matrix
%    are assumed zeros. This is s form of sparse matrix storage, useful if
%    the matrix consists predominantly of zeros.
% 
%    dims are the matrix' dimensions [ny], [ny nx], [ny, nx, nz] etc
%
% TO 100925;

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

X=zeros(dims);
X(IV(:,1))=IV(:,2);
