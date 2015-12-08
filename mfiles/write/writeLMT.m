function writeLMT(basename,lmt)
%WRITELMT writes input file for MT3MDS LMT package
%
% Example
%    writeLMT(basename,lmt) --- write LMT file
%
% TO 100101

fid=fopen([basename,'.LMT'],'wt');

%0.
fprintf(fid,'%s\n',['# MATLAB writeLMT ' datestr(now)]);
fprintf(    '%s\n',['# MATLAB writeLMT ' datestr(now)]);

fprintf(fid,'OUTPUT_FILE_NAME        %s\n',lmt.OUTPUT_FILE_NAME);
fprintf(fid,'OUTPUT_FILE_UNIT        %d\n',333);   % 333 is the default unit number for FTL
fprintf(fid,'OUTPUT_FILE_HEADER      %s\n',lmt.OUTPUT_FILE_HEADER);
fprintf(fid,'OUTPUT_FILE_FORMAT      %s\n',lmt.OUTPUT_FILE_FORMAT);

fclose(fid);
