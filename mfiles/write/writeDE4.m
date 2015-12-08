function writeDE4(basename,de4)
%WRITEDE4 writes input file for MODFLOW's direct sovler package DE4
%
% Example:
%    writeDE4(basename,de4) --- write DE4 file
%
% TO 100101

fid=fopen([basename,'.',de4.ext],'wt');

%0.
fprintf(fid,'%s\n',['# MATLAB writeDE4 ' datestr(now)]);
fprintf(    '%s\n',['# MATLAB writeDE4 ' datestr(now)]);

%F1 ITMX	MXUP	MXLOW	MXBW
%    ITMX max # of outer iteration, use >1
%    MXUP  max number of upper equations (use 0)
%    MXLOQ max number of lower equations (use 0)
%    MXBW  maximum bandwidth (use 0);
fprintf(fid,'%10d%10d%10d%10d     ITMX MXUP MXLOW MXBW\n',...
    de4.ITMX,de4.MXUP,de4.MXLOW,de4.MXBW);

%F2 IFREQ	MUTD4	ACCL  HCLOSE IPRD4
%   IFREQ  flag indicating at which coefficients in [A] change.
%      1 = flow equations are linear, coefs const for all stress periods
%      2 = flow equations are linear, but coefficients of [A] for some
%          stress periods may change at start of stress period.
%      3 = nonlnear flow equations -> some coefficients of [A] depend on
%          head (wetting etc). This requires ITMX>1
%   MUT4D flag indicating the quantity of information printed when printing
%         at each step.
%      0 = # of iterations in time step and maximum head change per iteration is printed
%      1 = only the # of iterations is printed
%      2 = no information is printed
%   ACCL is a multiplier for the computed head change for eah iteration.
%        Normally 1. Someimes a value >1 may be help convergence when using
%        external iteration (IFREQ=3) to solve nonlinear problems.
%        Always use 1 for linear problems.
%   HCLOSE head closure criterion if ITMX>1. Only used if nonliear.
%   IPRD4 time step interval for printing of convergence info when ITMX>1.
fprintf(fid,'%10d%10d%10.4g%10.4g%10d    IFREQ MUTD4 ACCL HCLOSE IPRD4\n',...
    de4.IFREQ,de4.MUTD4,de4.ACCL,de4.HCLOSE,de4.IPRD4);

fclose(fid);
