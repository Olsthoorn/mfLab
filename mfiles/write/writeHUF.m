function  writeHUF(basename,huf)
%WRITEHUF writes input file for MODFLOW's HUF package
%
% Example:
%    writeLPF(basename,lpf);
%
% ToDo: Untested, lots to be one here (090812)
%
% TO 090812

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fid=fopen([basename,'.',huf.ext],'wt');

%0.
fprintf(fid,'%s\n',['# MATLAB writeHUF ' datestr(now)]);
fprintf(    '%s\n',['# MATLAB writeHUF ' datestr(now)]);

fid=fopen([basename,'.',lpf.ext],'wt');

%1.
fprintf(fid,'%10d%10.3g%10d%10d%10d     IHUFCB HDRY NHUF NPHUF IOHUF\n',...
    huf.IHUFCB,huf.HDRY,huf.NPHUF,huf.IOHUF);


%2 Layer convertibility (confined / unconfined)
warray(fid,huf.LTHUF,huf.unit,'(25I3)','LTHUF',false);

%3 Layer wet-dry switch 
warray(fid,huf.LWYWT,huf.unit,'(25I3)','LAYWT',false);

%4  Defaults, don't touch, if at least one layer is wettable
WETFCT=1; IWETIT=1; IHDWET=0;
if any(huf.LAYWET)
    fprintf(fid,'%10.3g%10d%10d     WETFCT IWETIT IHDWET\n',WETFCT,IWETIT,IHDWET);
end

%5 WETDRY combination of layer threshold and flag indicating which
%neighboring cells can wet this cell
huf.WETDRY=zeros(huf.NROW,huf.NCOL,sum(huf.LAYWT>0));
k=0;
if any(huf.LAYWT)
    for i=1:huf.NLAY
        if huf.NLAY
            k=k+1;
            warray(fid,huf.WETDRY(iL)*PLANE,huf.unit,'(25I3)',sprintf(' WETDRY{%d}',iL),true,huf.FREE);  %Sy
        end
    end
end

%6
huf.TOP=zeros(huf.NROW,huf.NCOL,huf.NHUF);
for iHuf=1:huf.NHUF
    %6
    fprintf(fid,'%s\n',huf.HGUNAM{iHuf});
    %7
    warray(fid,huf.TOP(:,:,iHuf),huf.unit,'(10E12.3)',sprintf('TOP{%d}',iHuf),true,huf.FREE);
    %8
    warray(fid,huf.THCK(:,:,iHuf),huf.unit,'(10E12.3)',sprintf('THCK{%d}',iHuf),true,huf.FREE);
end

%9
for iHuf=1:huf.NHUF
    fprintf(fid,'%10s%10.4f%10.4f\n',HGUNAM,HGUHANI,HGUVANI);
    if strcmpi(HGUNAM,'ALL')
        for jHuf=1:huf.NHUF
            huf.HGUHANI(jHuf)=HGUHANI;
            huf.HGUVANI(jHif)=GHUVANI;
        end
        break;
    else
        found=0;
        for jHuf=1:huf.NHUF
            if strcmpi(HGUNAM,huf.HGUNAM(jHuf))
                huf.HGUHANI=HGUHANI;
                huf.HGUVANI=HGUVANI;
                found=1;
            end
        end
        if ~found
            error('HGUNAM=%s must be defined before reading HGUHANI and HGUVANI',HGUNAM);
        end
    end
end    

for iPar=1:huf.NPHUF
%10 PARNAM PARTYP Parval NCLU
    fprintf(fid,'%10s%10s%12.4f%10d\n',huf.PARNAM{iPar},huf.PARTYP{iPar},huf.Parval(iPar),huf.NCLU(iPar));
%11  parameter clusters
    for iClu=1:huf.NCLU(iPar)
        fprintf(fid,'%10s%10s%10s\n',huf.Cluster(iPar).HGUNAM{iClu},...
                                     huf.Cluster(iPar).Mltarr{iClu},...
                                     huf.Cluster(iPar).Zonarr{iClu});
        fprintf(fid,' %d',huf.Cluster(iPar).IZ);
        fprintf(fid,'\n');
    end
end

%12, optional printcodes and printflags
for iHuf=1:huf.NHUF
    fprintf(fid,'%10s%10s%10s%10s\n','PRINT',huf.HGUNAM{iHuf},huf.PRINTCODE{iHuf});
    fprintf(fid,'%10s',huf.PRINTFLAGS{iHuf});
    fprintf(fid,'\n');
end

fclose(fid);
