function writeSOR(basename,sor)
%WRITESOR writes input file for MODFLOW's SOR solver package.
%
% Example:
%    writeSOR(basename,sor) --- write SOR file
%
% TO 090101

fid=fopen([basename,'.',sor.ext],'wt');

%0.
fprintf(fid,'%s\n',['# MATLAB writeSOR ' datestr(now)]);
fprintf(    '%s\n',['# MATLAB writeSOR ' datestr(now)]);

%F1 MXITER
%    MXITER max # of iterations allowed in a time step
fprintf(fid,'%10d     MXITER\n',sor.MXITER);

%F2 ACCL  HCLOSE IPRSOR
%   ACCL accelleration variable usually between 1.0 and 2.0
%   HCLOSE it the head change criterion for convergence
%   IPRSOR time step print out interval voor SOR. If it is equal to zero
%          it is automatically set to 999.
fprintf(fid,'%10.4g%10.4g%10d    ACCL HCLOSE IPRSOR\n',...
    sor.ACCL,sor.HCLOSE,sor.IPRSOR);

fclose(fid);
