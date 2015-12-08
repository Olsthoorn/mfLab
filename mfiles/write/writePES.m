function writePES(basename,pes)
%WRITEPES writes iput file for PEST program (calibration for any model)
%
% Example:
%    writePES(basename,pes) ---- write pest file
%
% TO 100101

fid=fopen([basename,'.pes'],'wt');

%0.
fprintf(fid,'%s\r\n',['# MATLAB writePES ' datestr(now)]);
fprintf(    '%s\n'  ,['# MATLAB writePES ' datestr(now)]);

%1.
fprintf(fid,'%10g%10g%10g%10g\r\n',pes.MAXITER,pes.MAXCHANGE,pes.TOL,pes.SOSC);

%2.
fprintf(fid,'%10g%10g%10g%10g%10g%10g%10g%10g%10g\r\n',pes.IBEFLG, pes.IYCFLG,pes.IOSTAR,pes.NOPT,pes.NFIT,pes.SOSR,pes.RMAR,pes.RMARM,pes.IAP);

%3
fprintf(fid,'%10g%10g%10g\r\n',pes.IPRCOV,pes.IPRINT,pes.LPRINT);

%4
fprintf(fid,'%10g%10g%10g\r\n',pes.CSA, pes.FCONV,pes.LASTX);

%5
fprintf(fid,'%10g%10g%10g\r\n',pes.NPNG,pes.IPR,pes.MPR);

fclose(fid);
