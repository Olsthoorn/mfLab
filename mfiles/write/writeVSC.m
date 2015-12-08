function writeVSC(basename,vsc)
%WRITEVSC writes intput file for SEAWAT's vicosity flow package file
%
% Example:
%    writeVSC(basename,vsc);
%
% See also: write VDF
%
% TO 091006

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

spaces='          ';

fid=fopen([basename,'.',vsc.ext],'wt');

%0.
fprintf(fid,'%s\n',['# MATLAB writeVSC ' datestr(now)]);
fprintf(    '%s\n',['# MATLAB writeVSC ' datestr(now)]);

%1. MT3DMUFLG
fprintf(fid,'%10d%s%s%s%s     MT3DMUFLG\n',vsc.MT3DMUFLG,spaces,spaces,spaces,spaces);

%2 VSCMIN VSCMAX
% minimum and maximum fluid viscosity; use 0 if not limited by value
fprintf(fid,'%10g%10g%s%s%s     VISCMIN VISCMAX\n',...
    vsc.VISCMIN,vsc.VISCMAX,spaces,spaces,spaces);

% MT3DMUFLG ==0  --> read in VSC using %6 and %7
% MT3DMUFLG  >0  --> use MT3DMUFLG as species nr to compute viscosity
% MT3DMUFLG ==-1 --> compose viscosity using MT3D species for conc

if vsc.MT3DMUFLG>=0
    %3 MTMUSPEC DMUDC(1) CMUREF(1)
    if all(vsc.MTMUSPEC ~= vsc.MT3DMUFLG)
        error('Worksheet VSC: if MTMUSPEC>=0 one of the MTMUSPEC must match MT3DMUFLG (even if 0)');
    end
    i=find(vsc.MTMUSPEC==vsc.MT3DMUFLG,1,'first');

    fprintf(fid,'%10g%10g%10g%s%s     VISCREF DMUDC(1) CMUREF(1)\n',...
        vsc.VISCREF,vsc.DMUDC(i),vsc.CMUREF(i),spaces,spaces);

elseif vsc.MT3DMUFLG==-1
    %3a VISCREF
     fprintf(fid,'%10g%s%s%s%s     VSCREF\n',vsc.VISCREF,spaces,spaces,spaces,spaces);

    %3b  NSMUEOS MUTEMPOPT
    fprintf(fid,'%10d%10d%s%s%s     NSMUEOS,MUTEMPOPT\n',...
        vsc.NSMUEOS,vsc.MUTEMPOPT,spaces,spaces,spaces);
    for i=1:vsc.NSMUEOS
         %3c  MTMUSPC(NSMUEOS) DMUDC(NSMUEOS) CMUREF(NSMUEOS)
         fprintf(fid,'%10d%10g%10g%s%s     MTMUSPEC(%d),DMUDC(%d),CMUREF(%d)\n',...
             vsc.MTMUSPEC(i),vsc.DMUDC(i),vsc.CMUREF(i),spaces,spaces,i,i,i);
    end

    if vsc.MUTEMPOPT>0
            fprintf(fid,'%d ' ,vsc.MTMUTEMPSPEC);
            fprintf(fid,'%10g',vsc.AMUCOEF);
            fprintf(fid,'     MTMUTEMPSPEC  AMUCOEF(1:%d)\n',length(vsc.AMUCOEF));
    end
end

if vsc.MT3DMUFLG==0
    for iPer=1:vsc.NPER
        %4 INVISC
        fprintf(fid,'%10d%s%s%s%s     INVISC(%d)\n',...
            vsc.INVISC(iPer),spaces,spaces,spaces,spaces,iPer);

        if vsc.INVISC(iPer)>0
        %5 VISC(NCOL,NROW) -- U2DREL
            for iLay=1:vsc.NLAY
                warray(fid,vsc.VISC{iPer}(:,:,iLay),'(10G12.5)',sprintf('VSC(layer=%d)',iLay),true,vsc.FREE);
            end
        end
    end
end

fclose(fid);
