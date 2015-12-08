function writeCOC(basename,coc)
%WRITECOC writes intput for output control pacakge of CFP (Conduit FLow Package) version of MODFLOW
%
% Example:
%    writeCOC(basename,coc);
%
% TO 090715 090715 110513

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fprintf('# MATLAB writeCOC %s\n',datestr(now));
fid=fopen([basename,'.',coc.ext],'wt');

%0-2 COC total nodes used in simulation
fprintf(    '# MATLAB writeCFP %s\n',datestr(now));
fprintf(    '# %s\n','conduit output file');
fprintf(fid,'%s\n%10d\n','# total number of nodes',size(coc.NODE_NUMBERS,1));

%3-4 COC the number that indicated to each node
fprintf(fid,'%s\n','# item 3: number of each node ');
fprintf(fid,'%10d\n',coc.NODE_NUMBERS);

%5-6 COC number of time step for node head and node flow
fprintf(fid,'%s\n','# item 5: Time step for each node head and flow ');
fprintf(fid,'%10d\n',coc.N_NTS);

%7-8 COC total number of pipes
fprintf(fid,'%s\n','# item 7: total number of pipes ');
fprintf(fid,'%10d\n',size(coc.PIPE_NUMBERS,1));

%9-10 COC the nomber indicated to each pipe 
fprintf(fid,'%s\n','# item 9: number of each pipe ');
fprintf(fid,'%10d\n',coc.PIPE_NUMBERS);

%11-12 COC number of time step for pipe flow rate and Reynolds numbers 
fprintf(fid,'%s\n','# item 11: Time step for each pipe flow and Re number');
fprintf(fid,'%10d\n',coc.T_NTS);

fclose(fid);
