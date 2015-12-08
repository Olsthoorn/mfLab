function found=fileGrep(filename,expr)
%FILEGREP find expr in file using regexp(fileContents,expr)
%
% USAGE:
%    found = fileGrep(filename,expr)
%
% TO130312

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fid = fopen(filename,'r');

if fid<0, error('%s: Can''t open file <<%s>>\n',mfilename,filename); end

found = false;

while ~feof(fid)
    s = fgets(fid);
    fprintf('%s',s);
    if regexp(s,expr),
        found = true;
        break;
    end
end
