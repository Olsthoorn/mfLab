function  writeGCG(basename,gcg)
%WRITEGCG writes inptu for MT3DMS's GCG generalized conjugate gradient sovler package.
%
% Example:
%    writeGCG(basename,gcg)
%
% TO 0706030 081227

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fid=fopen([basename,'.',gcg.ext],'wt');

%C1 HEADING 1+2 (<=80 chars) -- No header allowed in GCG file
fprintf(    '%s\n',['# MT3DMS writeGCG ' datestr(now)]);


%F1 MXITER ITER1 ISOLVE NCRS (free)
%    MXITER max # of outer iteration, use >1
%    ITER1  max number of inner iterations (use 30-50)
%    ISOLVE preconditioner type
%       1= Jacobi
%       2= SSOR
%       3= modified incomplete choleski (MIC)
%    NCRS dispesion tensor treatment flag
%       0=lump all dispersion cross terms to RHS
%       1=include full dispersion tensor (memory intensive)
fprintf(fid,'%10d%10d%10d%10d     MXITER ITER1 ISOLVE NCRS\n',...
    gcg.MXITER,gcg.ITER1,gcg.ISOLVE,gcg.NCRS);

%F2 ACCL CCLOSE IPRGCG (free)
%   ACCL relaxtion fctor for SSOR optoin, use 1.0 as default
%   CCLOSE relative concentration convergence criterion e.g. 10e-4 to 10e-6
%   IPRGCG interval for printing maximu conc changes of each iteration
%     0=only print at end of stress period
fprintf(fid,'%10g%10g%10d               ACCL CCLOSE IPRGCG\n',...
    gcg.ACCL,gcg.CCLOSE,gcg.IPRGCG);

fclose(fid);
