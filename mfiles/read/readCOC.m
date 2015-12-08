function coc=readCOC(fname,pth)
%READCOC reads output control for CFP package (mf05)
%
% Example:
%    coc=readCOC([pth fname])
%
% TO 090715

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fprintf('# MATLAB readCOC %s\n',datestr(now));
fid=fopen([pth fname],'r');

% NODES ==========
fprintf(fgets(fid));                        % 0 required comment line
fprintf(fgets(fid));                        % 1 required comment line
coc.NNODES=fscanf(fid,'%d',1); fgets(fid);      % 2 NNODES
fprintf(fgets(fid));                        % 3 required comment line
coc.Nodes=zeros(coc.NNODES,1);                  % 4 Nodes
for i=1:coc.NNODES
    coc.Nodes(i)=fscanf(fid,'%d',1); fgets(fid);
end
fprintf(fgets(fid));                        % 5 required comment line
coc.N_NTS=fscanf(fid,'%d',1); fgets(fid);   % 6 time step interval for node head and flow output

% PIPES ===========
fprintf(fgets(fid));                        % 7 required comment line
coc.NPIPES=fscanf(fid,'%d',1); fgets(fid);      % 8 number of pipes
fprintf(fgets(fid));                        % 9 required comment line
coc.Pipes=zeros(coc.NPIPES,1);                  % 10 pipes
for i=1:coc.NPIPES
    coc.Pipes(i)=fscanf(fid,'%d',1); fgets(fid);
end
fprintf(fgets(fid));                    % 11 Required comment line
coc.T_NTS=fscanf(fid,'%d',1); fgets(fid);   % 12 pipe output time interval

fclose(fid);
