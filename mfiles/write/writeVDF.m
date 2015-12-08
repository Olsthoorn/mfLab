function writeVDF(basename,vdf)
%WRITEVDF writes input file for the variable density flow package file (SEAWAT)

% Example:
%    writeVDF(basename,vdf)
%
% TO 100101

fid=fopen([basename,'.',vdf.ext],'wt');

%0.
fprintf(fid,'%s\n',['# MATLAB writeVDF ' datestr(now)]);
fprintf(    '%s\n',['# SEAWAT writeVDF ' datestr(now)]);

%1. MTDNCONC MFADVFD NSWCPL IWTABLE
% MTDDRHOFLAG = MT3D species number that will be used
%     0=fluid density specified using input records 6 and 7
%     >0 fluid densisty computed using MT3D species number = MTDNCONC 
% MFADVFD  = flag internodal density values calculation method
%    2=central in space algoritm
%    <>2=upstream weighted algoritm
% NSWCPL = max number of non-linear coupling iterations
%  0 or 1  flow and transprot explicitly coupled using 1 one timestep lag
%  >1 implicit coupling using this number of coupling iterations
% IWTABLE = variable density water table correction computation flag
%    0=no correction
%    >0=correction is applied

spaces='          ';
fprintf(fid,'%10d%10d%10d%10d     MTD3DRHOFLG MFADVFD NSWTCPL IWTABLE\n',...
    vdf.MT3DRHOFLG,vdf.MFNADVFD,vdf.NSWTCPL,vdf.IWTABLE);

%2 DENSEMIN DENSEMAX
% minimum and maximum fluid densities use 0 if not limited by value
fprintf(fid,'%10g%10g%s%s     DENSEMIN DENSEMAX\n',...
    vdf.DENSEMIN,vdf.DENSEMAX,spaces,spaces);

if vdf.NSWTCPL>1 || vdf.NSWTCPL==-1
    %3  DNSCRIT convergence criterion   /// not in SW V4 !!!
     fprintf(fid,'%10g%s%s%s     DENSCRIT\n',vdf.DNSCRIT,spaces,spaces,spaces);
end

%4  DENSEREF DENSESLP
%    DENSEREF reference density of freshwater
%    DENSESLP density slopw drho/dc=20000/25 mgCl/L per 25 kg/m3
if vdf.MT3DRHOFLG>=0
    fprintf(fid,'%10g%10f%s%s     DENSEREF DRHODC\n',vdf.DENSEREF,vdf.DRHODC,spaces,spaces);
elseif vdf.MT3DRHOFLG==-1
    %4a
     fprintf(fid,'%10g%10g%12f%s     DENSEREF DRHODPRHD PRHDREF\n',vdf.DENSEREF,vdf.DRHODPRHD,vdf.PRHDREF,spaces);

    %4b  NSRHOEOS
    fprintf(fid,'%10d%s%s%s     NSRHOEOS\n',vdf.NSRHOEOS,spaces,spaces,spaces);
    for i=1:vdf.NSRHOEOS

         %3d  MTRHOSPEC(NSRHOEOS) DRHODC(NSRHOEOS) CRHOREF(NSRHOEOS)
         fprintf(fid,'%10d%10f%10g%s     MTRHOSPEC(%d),DRHODC(%d),CRHOREF(%d)\n',...
             vdf.MTRHOSPEC(i),vdf.DRHODC(i),vdf.CRHOREF(i),spaces,i,i,i);
    end
end
   
%5  FRSTDT first timestep length whenever IMT process is active
fprintf(fid,'%10g%s%s%s     FIRSTDT\n',vdf.FIRSTDT,spaces,spaces,spaces);

if vdf.MT3DRHOFLG==0
    for iPer=1:vdf.NPER
    %6 INDENSE
        fprintf(fid,'%10d\n',vdf.INDENSE(iPer));
        if vdf.INDENSE(iPer)>0
        %7 DENSE(NCOL,NROW) -- U2DREL
            for iLay=1:vdf.NLAY
                warray(fid,vdf.DENSE{iPer}(:,:,iLay),'(10G12.5)',sprintf('DENSE(layer=%d)',iLay),true,vdf.FREE);
            end

        end
    end
end


fclose(fid);
