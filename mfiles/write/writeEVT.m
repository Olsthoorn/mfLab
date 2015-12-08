function writeEVT(basename,evt,evtp)
%WRITEEVT writes input file for MODFLOW's EVT package
%
% Example:
%    writeEVT(fname,evt)  ---- write EVT file
%
% TO 091007

% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fid=fopen([basename,'.EVT'],'wt');

%0.
fprintf(fid,'# MATLAB writeEVT, %s\n',datestr(now));
fprintf(    '# MATLAB writeEVT, %s\n',datestr(now));

%1.
%if evtp.NPEVT~=0
    fprintf(fid,'PARAMETER%10d\n',[evtp.NPEVT]);
%end

%2.
fprintf(fid,'%10d%10d  NEVTOP  IEVTCB\n',[evt.NEVTOP,evt.IEVTCB]);

%3
if evtp.NPEVT~=0
    for nP=1:evtp.NPEVT
        fprintf(fid,'%11s%11s% 10d%10d\n',char(evtp.PARNAM(nP)),char(evtp.PARTYP(nP)),evtp.Parval(nP),evtp.NCLU(nP));
        for nC=1:evtp.NCLU(nP)
        %4b.
        fprintf(fid,'%11s%11s%10i\n',char(evtp.Mltarr(nP)),char(evtp.Zonarr(nP)),evtp.IZ(nP));
        end
    end
end

for iP=1:evt.NPER
    %5
    fprintf(fid,'%10d%10d%10d%10d     INSURF(%d) INEVTR(%d) INEXDP(%d) INIEVT(%d)\n',...
        evt.INSURF(iP),evt.INEVTR(iP),evt.INEXDP(iP),evt.INIEVT(iP),iP,iP,iP,iP);
    
    %6
    if evt.INSURF(iP)>=0
        if iscell(evt.SURF)
            warray(fid,evt.SURF{iP},evt.unit,'(10E15.6)',sprintf(' SURF{%d}',iP),true,evt.FREE);
        else
            warray(fid,evt.SURF(:,:,iP),evt.unit,'(10F15.6)',sprintf(' SURF{%d}',iP),true,evt.FREE);
        end
    end
    
    %7
    if evt.INEVTR(iP)>=0
        if iscell(evt.EVTR)
            warray(fid,evt.EVTR{iP},evt.unit,'(10F15.6)',sprintf(' EVTR{%d}',iP),true,evt.FREE);
        else
            warray(fid,evt.EVTR(:,:,iP),evt.unit,'(10F15.6)',sprintf(' EVTR{%d}',iP),true,evt.FREE);
        end
    end

    %9
    if evt.INEXDP(iP)>=0
        if iscell(evt.EXDP)
            warray(fid,evt.EXDP{iP},evt.unit,'(10F15.6)',sprintf(' EXDP{%d}',iP),true,evt.FREE);
        else
            warray(fid,evt.EXDP(:,:,iP),evt.unit,'(10F15.6)',sprintf(' EXDP{%d}',iP),true,evt.FREE);
        end
    end
    
    %10
    if evt.NEVTOP==2 && (evt.INIEVT(iP)>=0 || iP==1) % if NEVTOP==1, INIEVT is obligaory in first stress period
        if iscell(evt.IEVT)
            warray(fid,evt.IEVT{iP},evt.unit,'(20I5)',sprintf(' IEVT{%d}',iP),true,evt.FREE);
        else
            warray(fid,evt.IEVT(:,:,iP),evt.unit,'(20I5)',sprintf(' IEVT{%d}',iP),true,evt.FREE);
        end
    end

end
fclose(fid);
