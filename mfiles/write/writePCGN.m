function writePCGN(basename,pcgn)
%WRITEPCG writes input for MODLFOW's PCG solver package
%
% Example:
%    writePCGN(basename,pcgn) --- write PCGN file
%
% TO 151115

fid=fopen([basename,'.',pcgn.ext],'wt');

%0.
fprintf(fid,'%s\n',['# MATLAB writePCGN ' datestr(now)]);
fprintf(    '%s\n',['# MATLAB writePCGN ' datestr(now)]);

%Set1:  ITER_MO, ITER_MI, CLOSE_R, CLOSE_H (free)
fprintf(fid,' %9d %9d %9g %9g     ITER_MO, ITER_MI, CLOSE_R, CLOSE_H\n',...
    pcgn.ITER_MO, pcgn.ITER_MI, pcgn.CLOSE_R, pcgn.CLOSE_H);

%Set2: RELAX, IFILL, UNIT_PC, UNIT_TS
fprintf(fid,' %9g %9d %9d% 9d    RELAX, IFILL, UNIT_PC, UNIT_TS\n',...
    pcgn.RELAX, pcgn.IFILL, pcgn.UNIT_PC, pcgn.UNIT_TS);

%Set3: ADAMP, DAMP, DAMP_LB, RATE_D, CHGLIMIT
fprintf(fid,' %9g %9g %9g %9g %9g     ADAMP, DAMP, DAMP_LB, RATE_D, CHGLIMIT\n',...
    pcgn.ADAMP, pcgn.DAMP, pcgn.DAMP_LB, pcgn.RATE_D, pcgn.CHGLIMIT);

%Set4: ACNVG, CNVG_LB, MCNVG, RATE_C, IPUNIT
fprintf(fid,' %9g %9d %9d %9g %9d     ACNVG, CNVG_LB, MCNVG, RATE_C, IPUNIT\n',...
    pcgn.ACNVG, pcgn.CNVG_LB, pcgn.MCNVG, pcgn.RATE_C, pcgn.IPUNIT);

fclose(fid);
