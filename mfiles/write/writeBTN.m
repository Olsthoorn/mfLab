function  writeBTN(basename,btn)
%WRITEBTN writes input file for MT3MDS's basic transport package.
%
% Example:
%    writeBTN(basename,btn) --- write basic transport package file
%
% TO 0706030 081227 100827

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

%btn.unit=100; % see doc MT3DMS p97 incompatibility with UD2REL

fid=fopen([basename,'.',btn.ext],'wt');

ctrlRecDesired = true;

%A1 HEADING 1+2 (<=80 chars)
 fprintf(fid,'%s\n',['# MATLAB writeBTN ' datestr(now)]);
 fprintf(    '%s\n',['# MT3DMS writeBTN ' datestr(now)]);

%A2
 fprintf(fid,'%s\n','# Input file for BTN package used by MT3DMS and SEAWAT');

 
%A3. logical flags for major transport and solution options (dummy in recent MT3DMS versions)
fprintf(fid,'%10d%10d%10d%10d%10d%10d     NLAY NROW NCOL NPER NCOMP MCOMP\n',...
    btn.NLAY,btn.NROW,btn.NCOL,btn.NPER,btn.NCOMP,btn.MCOMP);

%A4 UNITS
 fprintf(fid,'%4s%4s%4s     TUNIT LUNIT MUNIT\n',btn.TUNIT,btn.LUNIT,btn.MUNIT);

%%5 TRNOP transport packages in use (flags T or F)
for i=1:numel(btn.TRNOP)
    if btn.TRNOP(i), fprintf(fid,' T'); else fprintf(fid,' F'); end
end
fprintf(fid,'     TRNOP(10): (ADV DSP SSM RCT GCG STR XXX XXX XXX XXX\n');

%%A6 LAYCON(NLAY) (40I2)
warray(fid,btn.LAYCON(:)',btn.unit,'(40I2)','LAYCON',~ctrlRecDesired);

%A7 DELC(NCOL) RARRAY reader (Note that conversion to rowvector is a
%must for correct printing!')
warray(fid,btn.DELR(:)',btn.unit,'(10E12.3)','DELR');

%A8 DELC(NROW) RARRAY reader (Not conversion to row vector for correct printing
warray(fid,btn.DELC(:)',btn.unit,'(10E12.3)','DELC');

%A9 HTOP(NCOL,NROW), RARRAY reader
warray(fid,btn.Z(:,:,1),btn.unit,'(10E12.3)','ZTOP');

%%A10 DZ(NCOL,NROW) RARRAY
for i=1:btn.NLAY
    warray(fid,btn.Z(:,:,i)-btn.Z(:,:,i+1),btn.unit,'(10E12.3)',sprintf('DZ{%d}',i));
end

%%A11 PRSITY(NCOL,NROW) RARRAY effective porosity (single porsoity system)
for i=1:btn.NLAY
    if numel(btn.PRSITY)==1
        warray(fid,btn.PRSITY       ,btn.unit,'(10E12.3)',sprintf('PRSITY{%d}',i));
    else
        warray(fid,btn.PRSITY(:,:,i),btn.unit,'(10E12.3)',sprintf('PRSITY{%d}',i));
    end
end

%%A12 ICBUND(NCOL,NROW) IARRAY
for i=1:btn.NLAY
    warray(fid,btn.ICBUND(:,:,i),btn.unit,'(20I4)',sprintf('ICBUND{%d}',i));
end

%%A13 STCONC(NCOL,NROW) RARRAY (one for each layer) start concentrations
for iComp=1:btn.NCOMP
    for iLay=1:btn.NLAY
       if iscell(btn.STCONC)
           warray(fid,btn.STCONC{iComp}(:,:,iLay),btn.unit,'(10E12.3)',sprintf('STCONC{iComp=%d,iLay=%d}',iComp,iLay));
       elseif numel(btn.STCONC)==1 
           warray(fid,btn.STCONC          ,btn.unit,'(10E12.3)',sprintf('STCONC{iComp=%d,Lay=%d}',iComp,iLay));
       else
           warray(fid,btn.STCONC(:,:,iLay),btn.unit,'(10E12.3)',sprintf('STCONC{iComp=%d,Lay=%d}',iComp,iLay));
       end
    end
end

%%A14 CINACT, THKMIN 2F10.0
fprintf(fid,'%10.3g%10.3g     CINACT THKMIN\n',...
    btn.CINACT,btn.THKMIN);

%%A15 IFMTCN IFMTNP IFMTRF IFMTDP SAVUCN  4I10 L10
if btn.SAVUCN, btn.SAVUCN='T'; else btn.SAVUCN='F'; end
fprintf(fid,'%10d%10d%10d%10d%10s     IFMTCN IFMTNP IFMTRF IFMTDP SAVUCN\n',...
    btn.IFMTCN,btn.IFMTNP,btn.IFMTRF,btn.IFMTDP,btn.SAVUCN);

%%A16 NPRS I10 ouput frequency flag
if isfield(btn,'TIMPRS')
    if btn.NPRS>0 && btn.TIMPRS(1)==0
        btn.TIMPRS(btn.TIMPRS==0)=[];  % zeros are not allowed!
        btn.NPRS  =numel(btn.TIMPRS);
    end
end
fprintf(fid,'%10d     NPRS\n',btn.NPRS);



%%A17 TIMPRS(NPRS) 8F10.0 total elapsed timeto print simulation resuls
if btn.NPRS>0
    for i=1:numel(btn.TIMPRS)  % squeeze (round) to make sure numbers are 10 wide
        btn.TIMPRS(i)=sscanf(numsqueeze(btn.TIMPRS(i),'%10g',10),'%g',1);
    end
    warray(fid,btn.TIMPRS,btn.unit,'(8G10.4)',sprintf('TIMPRS(1..%d)',numel(btn.TIMPRS)),'norec');
end

%%A18 NOBS, NPROBS 2I10
btn.NOBS=size(btn.OBS,1);
fprintf(fid,'%10d%10d     NOBS NPROBS\n',btn.NOBS,btn.NPROBS);
%%A19 KOBS IOBS JOBS 3I10
if btn.NOBS>0
    for iObs=1:btn.NOBS
        fprintf(fid,'%10d%10d%10d     L R C obs point(%d)\n',...
            btn.OBS(iObs,1:3),iObs);
    end
end

%%A20 CHKMAS NPRMAS  L10 I10
if btn.CHKMAS, btn.CHKMAS='T'; else btn.CHKMAS='F'; end
fprintf(fid,'%10s%10d     CHKMAS NPRMAS\n',btn.CHKMAS,btn.NPRMAS);

% for each stress period
for i=1:btn.NPER
    % Make sure we always geometrial increase of time step
    if btn.TSMULT(i)<=0, btn.TSMULT=1.25; end
    
    %%A21 PERLEN, NSTP, TSMULT F10.0 I10 F10.0
    fprintf(fid,'%10s%10d%10.3g               PERLEN(%d) NSTP(%d) TLMULT(%d)\n',...
        numsqueeze(btn.PERLEN(i),'%10g',10),btn.NSTP(i),btn.TSMULT(i),i,i,i);

    %%A22 if TSMULT<=0 --> TSLNGH(NSTP) -- skipped  only geometrical time steps
    % We force TSMULT to be >= 1 sie above !!
    
    %%A23 D10, MXSTRN, TTSMULT, TTSMAX F10.0 I10 2F10.0
    fprintf(fid,'%10.3g%10d%10.3g%10.3g     DT0(%d) MXSTRN(%d) TTSMULT(%d) TTSMAX(%d)\n',...
        btn.DT0(i),btn.MXSTRN(i),btn.TTSMULT(i),btn.TTSMAX(i),i,i,i,i);

end

fclose(fid);
