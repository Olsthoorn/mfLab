function writeETS(basename,ets,etsp)
%WRITEETS writes input file for MODFLOW's ETS package
%
% Example:
%    writeEVT(fname,ets)  ---- write ETS file
%
% TO 091007 120312

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fid=fopen([basename,'.ETS'],'wt');

%0.
fprintf(fid,'# MATLAB writeETS, %s\n',datestr(now));
fprintf(    '# MATLAB writeETS, %s\n',datestr(now));

%1.
fprintf(fid,'%10d%10d%10d%10d  NETSOP  IETSCB NPETS NETSEG\n',ets.NEVTOP,ets.IEVTCB,etsp.NPEVT,ets.NETSEG);

%2+%3
if etsp.NPEVT~=0
    for nP=1:etsp.NPEVT
        %2
        fprintf(fid,'%11s%11s% 10d%10d\n',char(etsp.PARNAM(nP)),char(etsp.PARTYP(nP)),etsp.Parval(nP),etsp.NCLU(nP));
        for nC=1:etsp.NCLU(nP)
        %3.
        fprintf(fid,'%11s%11s%10i\n',char(etsp.Mltarr(nP)),char(etsp.Zonarr(nP)),etsp.IZ(nP));
        end
    end
end

for iP=1:ets.NPER
    %4
    fprintf(fid,'%10d%10d%10d%10d%10d     INSURF(%d) INEVTR(%d) INEXDP(%d) INIEVT(%d) INSGDF(%d)\n',...
        ets.INSURF(iP),ets.INEVTR(iP),ets.INEXDP(iP),ets.INIEVT(iP),ets.INSGDF(iP),iP,iP,iP,iP,iP);
    
    %5 -- SURF == ETSS (surface elevation above which ET is max)
    if ets.INSURF(iP)>=0
        if iscell(ets.SURF)
            warray(fid,ets.SURF{iP},ets.unit,'(10E15.6)',sprintf(' SURF{%d}',iP),true,ets.FREE);
        else
            warray(fid,ets.SURF(:,:,iP),ets.unit,'(10E15.6)',sprintf(' SURF{%d}',iP),true,ets.FREE);
        end
    end
    
    %6 -- EVTR == ETSR (maximum evapotranspiration)
    if ets.INEVTR(iP)>=0
        if iscell(ets.EVTR)
            warray(fid,ets.EVTR{iP},ets.unit,'(10E15.6)',sprintf(' EVTR{%d}',iP),true,ets.FREE);
        else
            warray(fid,ets.EVTR(:,:,iP),ets.unit,'(10E15.6)',sprintf(' EVTR{%d}',iP),true,ets.FREE);
        end
    end

    %8 -- ETSX == EXDP (extension depth)
    if ets.INEXDP(iP)>=0
        if iscell(ets.EXDP)
            warray(fid,ets.EXDP{iP},ets.unit,'(10E15.6)',sprintf(' EXDP{%d}',iP),true,ets.FREE);
        else
            warray(fid,ets.EXDP(:,:,iP),ets.unit,'(10E15.6)',sprintf(' EXDP{%d}',iP),true,ets.FREE);
        end
    end
    
    %9 -- IETS == IEVT (layer from which ET extracts water)
    if ets.NEVTOP>1 && (ets.INIEVT(iP)>=0 || iP==1) % if NEVTOP==1, INIEVT is obligaory in first stress period
        if iscell(ets.IEVT)
            warray(fid,ets.IEVT{iP},ets.unit,'(20I5)',sprintf(' IEVT{%d}',iP),true,ets.FREE);
        else
            warray(fid,ets.IEVT(:,:,iP),ets.unit,'(20I5)',sprintf(' IEVT{%d}',iP),true,ets.FREE);
        end
    end
    
    %10+%11
    if ets.INSGDF(iP)>=0 % then segments are defined
        for iSeg=1:ets.NETSEG-1
            %10 -- fraction of extention depth where this segment ends)
            if iscell(ets.IEVT)
                warray(fid,ets.PXDP{iP},ets.unit,'(10E15.6)',sprintf(' PXDP{%d}',iP),true,ets.FREE);
            else
                warray(fid,ets.PXDP(:,:,iP),ets.unit,'(10E15.6)',sprintf(' PXDP{%d}',iP),true,ets.FREE);
            end
            %11 -- fraction of max ET for this point, i.e. section end
            if iscell(ets.PETM)
                warray(fid,ets.PETM{iP},ets.unit,'(10E15.6)',sprintf(' PETM{%d}',iP),true,ets.FREE);
            else
                warray(fid,ets.PETM(:,:,iP),ets.unit,'(10E15.6)',sprintf(' PETM{%d}',iP),true,ets.FREE);
            end
        end
    end
end
fclose(fid);
