function writePCG(basename,pcg)
%WRITEPCG writes input for MODLFOW's PCG solver package
%
% Example:
%    writePCG(basename,pcg) --- write PCG file
%
% TO 100101

fid=fopen([basename,'.',pcg.ext],'wt');

%0.
fprintf(fid,'%s\n',['# MATLAB writePCG ' datestr(now)]);
fprintf(    '%s\n',['# MATLAB writePCG ' datestr(now)]);

%F1 MXITER ITER1 NPCOND (free)
%    MXITER max # of outer iteration, use >1
%    ITER1  max number of inner iterations (use 30-50)
%    NPCOND conditioner
%       1= modified incomplete Cholesky
%       2= polynomica (use for vector computers)
fprintf(fid,'%10d%10d%10d%40s     MXITER ITER1 NPCOND\n',...
    pcg.MXITER,pcg.ITER1,pcg.NPCOND,' ');

%F2 HCLOSE RCLOSE RELAX NBPOL IPRPCG MUTPCG DAMP
%   HCLOSE head change criterion for convergence
%   RCLOSE residulal change criterion for convergence (L^3/T)
%   RELAX  Relaxation perameter used with NPCOND==1. Usually RELAX==1, but
%      for some porblems a value of 0.99, 0.98 or 0.97 improves convegence.
%   NBPOL  is used when NPCOND==2 to indicate whether the estimate of the
%        upper bound of the maximum eigen value is 2.0 of whether the estimate
%        will be calculated. NPBOL is used to specify the value is 2.0. For
%        any other value the estimate will be calculated.
%   CCLOSE relative concentration convergence criterion e.g. 10e-4 to 10e-6
%   IPRPCG printout interval for PCG. If set to 0 it will be reset to 999.
%   MUTPCG Flag controlling printing convergence information.
%      0 = print tables and max head change and residual each iteration
%      1 = print total number of iterations
%      2 = no printing
%      3 = print only if convergence fails
%   DAMP   Damping factor 0<DAMP<=1, typically DAMP=1, indicating damping.
fprintf(fid,'%10.4g%10.4g%10.4g%10d%10d%10d%10.4g    HCLOSE RCLOSE RELAX NBPOL IPRPCG MUTPCG DAMP\n',...
    pcg.HCLOSE,pcg.RCLOSE,pcg.RELAX,pcg.NBPOL,pcg.IPRPCG,pcg.MUTPCG,pcg.DAMP);

fclose(fid);
