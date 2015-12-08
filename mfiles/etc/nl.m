function nl(fid)
%NL print a new line
%
% USAGE:
%    nl(fid);
%
% TO 140415

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if nargin==0, fprintf('\n');
else          fprintf(fid,'\n');
end