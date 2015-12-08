function fprintfs(varargin)
%FPRINTFS fprintf but for multiple strings given in cell array
%
% USAGE:
%    fprintfs(fid,fmt,strs); 
%
%   fid and fmt as usual, but must contain %s to print a string.
%   strs a cell array of strings to print in sequence using a single %s.
%
%   Example:
%      fprintfs('<<%s>>',{'set','of','strings'})
%      fprintfs(' %s',{'set','of','strings'})
%      fprintfs(' Hello: %s,',{'set','of','strings'})
%      fprintfs(fid,'<<%s>>',{'set','of','strings'})
%
% TO 121124 130403

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if nargin>2
    [fid,varargin] = getNext(varargin,'double',1);
else
    fid = 1;
end
[fmt,varargin] = getNext(varargin,'char',' <<%s>>');
[strs, ~     ] = getNext(varargin,'cell',[]);
if isempty(strs)
    error('%s: last argument must be a cell array of strings');
end

for i=1:numel(strs)-1
    fprintf(fid,fmt,strs{i});
end
fprintf(fid,fmt,strs{end});
