function mult=readMULT(fname,mult)
%READMULT reads MODFLOW's multiplier file package
%
% Example:
%    readMUL(basename,mult);
%
% TO 070630 090713 090718

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fprintf('# MATLAB readMULT %s\n',datestr(now));

fid=fopen(fname,'r');

skipmodflowcomments(fid)

%1
s=fgets(fid);
mult.NArray=sscanf(s,'%d',1);

for i=1:mult.NArray
    %2
    s=fgets(fid); mult.name{i}=sscanf(s,'%s',1);
    %3
    mult.A{i}=mudread(fid,[mult.NROW,mult.NCOL]);
end

fclose(fid);
