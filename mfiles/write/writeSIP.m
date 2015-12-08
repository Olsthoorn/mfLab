function writeSIP(basename,sip)
%WRITESIP writes intput for MODFLOW's SIP solver package
%
% Example:
%    writeSIP(basename,sip)
%
% TO 100101

fid=fopen([basename,'.',sip.ext],'wt');

%0.
fprintf(fid,'%s\n',['# MATLAB writeSIP ' datestr(now)]);
fprintf(    '%s\n',['# MATLAB writeSIP ' datestr(now)]);

%F1 MXITER NPARM
%    MXITER max # of time through iteration loop in one time step.
%    NPARM the number of iteration variables to be used.
fprintf(fid,'%10d%10d     MXIter NPARM\n',sip.MXITER,sip.NPARM);

%F2 ACCL  HCLOSE IPCALC WSEED IPRSIP
%   ACCL accelleration variable > 0 and usually equal to 1.
%   HCLOSE the head change criterion for convergence.
%   IPCALC flag indicating where the seed for calculating iteration variables
%          will come from.
%     0 = the seed entered by user will be used.
%     1 = the seed will be calculated at the start os the simulation from
%         problem variables.
%   WSEED is the seed for calculating iteration variables. It is always
%   read, but it is only used if IPCALC==0.
%   IPRSIP printout interval for SIP. If 0 it will be set to 999 automatiically.
fprintf(fid,'%10.4g%10.4g%10d%10.4g%10d    ACCL HCLOSE IPCALC WSEED IPRSIP\n',...
    sip.ACCL,sip.HCLOSE,sip.IPCALC,sip.WSEED,sip.IPRSIP);

fclose(fid);
