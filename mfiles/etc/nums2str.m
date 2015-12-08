function s = nums2str(values,fmt1,fmt2)
%NUMS2STR converts numeric values to strings using given formats
%
% USAGE:
%    s = nums2str(values,fmt1[,fmt2]);
%
%    put vallues to string according to 
%    fmt1 (for the individual values) and fmt (to finish off, typically '\n'
%
% TO 121120

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fmt = repmat(fmt1,[1 numel(values)]);
if nargin<2, fmt2=''; end
    
s = sprintf([fmt fmt2],values);
    