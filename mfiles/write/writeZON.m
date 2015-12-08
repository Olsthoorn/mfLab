function writeZON(basename,zon)
%WRITEZON writes input for zone file
%
% Example:
%    writeZONE(basename,zone);
%
% TO 090101

fid=fopen([basename,'.zon'],'wt');

%0.
fprintf(fid,'%s\r\n',['# MATLAB writeZONE ' datestr(now)]);
fprintf(    '%s\n'  ,['# MATLAB writeZONE ' datestr(now)]);

%1.
fprintf(fid,'%10i\r\n',[zon.NZN]);
    
for iP=1:zon.NZN
    %2
    fprintf(fid,'%s\r\n',zon.ZONNAM{iP});
    
    %3
    warray(fid,zon.IZON{iP},'% g');
end
fclose(fid);

