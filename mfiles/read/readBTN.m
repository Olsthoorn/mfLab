function  btn=readBTN(basename,btn)
%READBTN reads MT3DMS/SEAWAT basic transport package file
%
% USAGE:
%    readBTN(basename,btn);
%
% TO 0706030 081227 100827

% Copyright 2007-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

btn.unit=100; % see doc MT3DMS p97 incompatibility with UD2REL

fid=fopen(basename,'rt');
if fid<0
    fid=fopen([basename,'.BTN'],'rt');
end

%A1 HEADING 1+2 (<=80 chars)
fprintf('%s\n',['# MATLAB readBTN ' datestr(now)]);
btn.Heading{1}=fgets(fid);
btn.Heading{2}=fgets(fid);
 
%A2
 fprintf('%s\n','# Input file for BTN package used by MT3DMS and SEAWAT');

 
%A3. logical flags for major transport and solution options (dummy in recent MT3DMS versions)
btn.NLAY  = fscanf(fid,'%10d',1); NLAY=btn.NLAY;
btn.NROW  = fscanf(fid,'%10d',1); NROW=btn.NROW;
btn.NCOL  = fscanf(fid,'%10d',1); NCOL=btn.NCOL;
btn.NPER  = fscanf(fid,'%10d',1);
btn.NCOMP = fscanf(fid,'%10d',1);
btn.MCOMP = fscanf(fid,'%10d',1);
fgets(fid);

%A4 UNITS
btn.TUNIT = fscanf(fid,'%4s',1);
btn.LUNIT = fscanf(fid,'%4s',1);
btn.MUNIT = fscanf(fid,'%4s',1);
fgets(fid);

%%5 TRNOP transport packages in use (flags T or F)
s=upper(fgetl(fid));
btn.TRNOP=sscanf(s,'%c'); btn.TRNOP(btn.TRNOP==' ')=''; %  TRNOP(10): (ADV DSP SSM RCT GCG STR XXX XXX XXX XXX\n');

%%A6 LAYCON(NLAY) (40I2)
btn.LAYCON=rarray(fid,NLAY,'norec');

%A7 DELC(NCOL) RARRAY reader (Note that conversion to rowvector is a
%must for correct printing!')
btn.DELR=rarray(fid,[1,NCOL]); % DELR

%A8 DELC(NROW) RARRAY reader (Not conversion to row vector for correct printing
btn.DELC=rarray(fid,[NROW,1]);  % DELC

%A9 HTOP(NCOL,NROW), RARRAY reader
btn.HTOP=rarray(fid,[NROW NCOL]);   % ZTOP

%%A10 DZ(NCOL,NROW) RARRAY
btn.DZ=NaN(NROW,NCOL,NLAY);
for iLay=1:btn.NLAY
    btn.DZ(:,:,iLay) = rarray(fid,[NROW, NCOL]); % DZ{%d}',i));
end

%%A11 PRSITY(NCOL,NROW) RARRAY effective porosity (single porsoity system)
btn.PRSITY=NaN(NROW,NCOL,NLAY);
for iLay=1:btn.NLAY
    btn.PRSITY(:,:,iLay)=rarray(fid,[NROW, NCOL]);  % PRSITY
end

%%A12 ICBUND(NCOL,NROW) IARRAY
btn.ICBUND=zeros(NROW,NCOL,NLAY);
for iLay=1:btn.NLAY
    btn.ICBUND(:,:,iLay)=rarray(fid,[NROW,NCOL]); %ICBUND
end

%%A13 STCONC(NCOL,NROW) RARRAY (one for each layer) start concentrations

btn.STCONC=cell(btn.NCOMP);
for iComp=1:btn.NCOMP
    for iLay=1:btn.NLAY
           btn.STCONC{iComp}(:,:,iLay)=rarray(fid,[NROW,NCOL]);  % STCONC
    end
end

%%A14 CINACT, THKMIN 2F10.0
btn.CINACT = fscanf(fid,'%f',1);
btn.THKMIN = fscanf(fid,'%f',1);
fgets(fid);

%%A15 IFMTCN IFMTNP IFMTRF IFMTDP SAVUCN  4I10 L10
%if btn.SAVUCN, btn.SAVUCN='T'; else btn.SAVUCN='F'; end
%fprintf('%10d%10d%10d%10d%10s     IFMTCN IFMTNP IFMTRF IFMTDP SAVUCN\n',...

s=fgetl(fid);
btn.IFMTCN= sscanf(s,'%d',1);
btn.IFMTNP= sscanf(s,'%d',1);
btn.IFMTRF= sscanf(s,'%d',1);
btn.IFMTDP= sscanf(s,'%d',1);
btn.SAVUCN= sscanf(s,'%s',1);

%%A16 NPRS I10 ouput frequency flag
s=fgetl(fid);
btn.NPRS=sscanf(s,'%d');
if btn.NPRS>0,
    btn.TIMPRS=rarray(fid,btn.NPRS,'norec')';
end % TIMPRS(1..NPRS)

%% A18 NOBS NPROBS
btn.NOBS   =fscanf(fid,'%d',1);
btn.NPROBS =fscanf(fid,'%d',1);
fgets(fid);

%% A19 KOBS IOBS JOBS
if btn.NOBS>0
    btn.OBS    = fscanf(fid,'%d',[3,btn.NOBS])'; 
end

%%A20 CHKMAS NPRMAS  L10 I10
btn.CHKMAS = fscanf(fid,'%s',1);
btn.NPRMAS = fscanf(fid,'%d',1);
fgets(fid);

% for each stress period
for i=1:btn.NPER
    
    %%A21 PERLEN, NSTP, TSMULT F10.0 I10 F10.0
    %fprintf('%10s%10d%10.3g 
    btn.PERLEN(i) = fscanf(fid,'%d',1);
    btn.NSTP(i)   = fscanf(fid,'%d',1); 
    btn.TSMULT(i) = fscanf(fid,'%f',1);
    fgets(fid);

    %%A22 length of time step within current time step if TMMULT<=0;
    if btn.TSMULT(i)<=0
        btn.TSLING{i}=fscanf(fid,'%f',btn.NSTP(i)); % time steps
        fgets(fid);
    end
            
    %%A23 D10, MXSTRN, TTSMULT, TTSMAX F10.0 I10 2F10.0
    %fprintf('%10.3g%10d%10.3g%10.3g     DT0(%d) MXSTRN(%d) TTSMULT(%d) TTSMAX(%d)\n',...
    btn.DT0(i)     = fscanf(fid,'%f',1);
    btn.MXSTRN(i)  = fscanf(fid,'%d',1);
    btn.TTSMULT(i) = fscanf(fid,'%f',1);
    btn.TTSMAX(i)  = fscanf(fid,'%f',1);
    fgets(fid);
end

fclose(fid);
