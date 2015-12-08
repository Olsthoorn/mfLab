function n = cols(X)
%COLS returns number of columns of array
%
% USAGE:
%    n=cols(X)
%
%    used with cellfun(@cols,...) in mf_setup (line 1829)
%
% TO 130411

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later
n = size(X,2);
