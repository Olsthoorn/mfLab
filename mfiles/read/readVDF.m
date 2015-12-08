function vdf=readVDF(fname,vdf)
%READVDF reads SEAWAT's VDF package
%
% Example:
%    writeVDF(fname,vdf)
%
% TO 090101

fid=fopen(fname,'r');

%0.
fprintf('# MATLAB writeVDF %s\n',datestr(now));

%1. MTDNCONC MFADVFD NSWCPL IWTABLE
% MTDNCONC = MT3D species number that will be used
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

vdf.MTDNCONC= fscanf(fid,'%f',1);
vdf.MFADVFD=  fscanf(fid,'%f',1);
vdf.NSWTCPL=  fscanf(fid,'%f',1);
vdf.IWTABLE=  fscanf(fid,'%f',1);
fgets(fid);

%2 DENSEMIN DENSEMAX
% minimum and maximum fluid densities use 0 if not limited by value
vdf.DENSEMIN=fscanf(fid,'%f',1);
vdf.DENSEMAX=fscanf(fid,'%f',1);
fprintf(fgets(fid));

% DENSCRIT
if vdf.NSWTCPL>1 || vdf.NSWTCPL==-1
    vdf.DENSCRIT=fscanf(fid,'%f',1);
    fgets(fid);
end

%4  DENSEREF DENSESLP
%    DENSEREF reference density of freshwater
%    DENSESLP density slopw drho/dc=20000/25 mgCl/L per 25 kg/m3
vdf.DENSEREF=fscanf(fid,'%f',1);
vdf.DENSESLP=fscnaf(fid,'%f',1);
fprintf(fgets(fid));

%5  FRSTDT first timestep length whenever IMT process is active
vdf.FIRSTDT=fscanf(fid,'%f',1);
fprintf(fgets(fid));

vdf.DENSE=cell(vdf.NPER,1);
vdf.CONC =cell(vdf.NPER,1);
for iP=1:vdf.NPER
    %6 INDENSE
    if vdf.MTDNCONC==0
        vdf.INDENSE(iP)=fscanf(fid,'%d',1);
        if vdf.INDENSE(iP)==1
            for iL=1:vdf.NLAY
                vdf.DENSE{iL}=mudread(fid,[vdf.NROW,vdf.NCOL]);
            end
        else
            for iL=1:vdf.NLAY
                vdf.CONC{iL}=mudread(fid,[vdf.NROW,vdf.NCOL]);
            end
        end
    end
end

fclose(fid);
