function s = sprintfs(fmt,strs)
%SPRINTFS sprintf for strings in cell array strs
%
% USAGE:
%    s = sprintfs(fmt,strs);
%
% TO 121124

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if ischar(strs)
    strs={strs};
end

s= char(zeros(1,255));
k=0;
for i=1:numel(strs)
    ss = sprintf(fmt,strs{i});
    s(k+(1:numel(ss))) = ss;
    k=k+numel(ss);
end
s=s(1:k);
