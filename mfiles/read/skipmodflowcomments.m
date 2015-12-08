function skipmodflowcomments(fid)
%SKIPMODFLOWCOMMENTS skips matlab comment lines, i.e. lines that start with #
%
% Example:
%    skipmodflowcomments(fid)
%    fid is a file indicator like cgf.
%
% TO 090713

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later
while 1
    p=ftell(fid);
    s=fgets(fid);
    if s(1)=='#',
        fprintf('%s',s);
    else
        fseek(fid,p,'bof');
        break;
    end
end
