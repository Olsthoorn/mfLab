function writeSEN(basename,sen)
%WRITESEN writes sensitivity file MODFLOW 2000
%
% Example:
%    writeSEN(basename,sen) ---- write sen file
%
%Tijs Dekker 070716

fid=fopen([basename,'.sen'],'wt');

%0.
fprintf(fid,'%s\r\n',['# MATLAB writeSEN ' datestr(now)]);
fprintf(    '%s\n'  ,['# MATLAB writeSEN ' datestr(now)]);

%1.
fprintf(fid,'%10g%10g%10g%10g\r\n',sen.NPLIST,sen.ISENALL,sen.IUHEAD,sen.MXSEN);

%2.
fprintf(fid,'%10g%10g%10g%10g\r\n',sen.IPRINTS, sen.ISENSU,sen.ISENPU,sen.ISENFM);

for iP=1:sen.NPLIST
%3
    fprintf(fid,'%10s%10g%10g%10g%10g%10g%10g\r\n',sen.PARNAM{iP},sen.ISENS(iP),sen.LN(iP),sen.B(iP),sen.BL(iP),sen.BU(iP),sen.BSCAL(iP));
end



fclose(fid);
