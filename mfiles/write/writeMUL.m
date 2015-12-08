function writeMUL(basename,mul)
%WRITEMUL writes input file for MODFLOW's multiplier array package
%
% Example:
%    writeMUL(basename,mul)  ---- write the multiplier file
%
% TO 100101

fid=fopen([basename,'.mul'],'wt');

%0.
fprintf(fid,'%s\r\n',['# MATLAB writeMUL ' datestr(now)]);
fprintf(    '%s\n'  ,['# MATLAB writeMUL ' datestr(now)]);

%1.
fprintf(fid,'%10i\r\n',[mul.NML]);
    
for iP=1:mul.NML
    %2
    fprintf(fid,'%10s\r\n',mul.MLTNAM{iP});
    
    %3
    warray(fid,mul.RMLT{iP});
end
fclose(fid);


