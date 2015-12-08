function  writeSSM(basename,ssm)
%WRITESSM writes input file for MT3DMS's sink-source mixing package SSM package
%
% Example:
%    writeSSM(basename,ssm)
%
% TO 0706030 081227

% Copyright 2009-2012 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

AsPrevious = -1;
Equal=zeros(ssm.NCOMP,1);

ssm.unit=100;  % incompatibility with UD2REL (MT3DMS p97)

fid=fopen([basename,'.',ssm.ext],'wt');

%D0 HEADING 1+2 (<=80 chars) -- No header allowd in DSP file
fprintf(    '%s\n',['# MT3DMS writeSSM ' datestr(now)]);

% Option flags for WEL, DRN etc packages + FNEW (n=1:4) future options
% See Zheng & Wang (1999), item D1, p119
FWEL=1;   % column nr of option flag 1
FDRN=2;   % column nr of option flag 2
FRCH=3;   % etc.
FEVT=4;
FRIV=5;
FGHB=6;
FCHD=7;
FMNW=8;

optionflag = zeros(1,10); % Ten flags defined in MT3DMS, initially false

optionflag (FRCH) =      ssm.FRCH;             % already determined form mf_setup
optionflag (FEVT) =      ssm.FEVT;             % same

if ~isfield(ssm,'PNTSRC'),
    ssm.PNTSRC = {};
else
    % column cell array
    ssm.PNTSRC = ssm.PNTSRC(:);
end

%% Add wells
if isfield(ssm,'well') && ~isempty(ssm.well)
    ssm.PNTSRC = [ssm.PNTSRC; ssm.well.PNTSRC()];
end

%% Add MNW1 wells
if isfield(ssm,'MNW1') && ~isempty(ssm.MNW1)
    ssm.PNTSRC = [ssm.PNTSRC; ssm.MNW1.PNTSRC()];
end
%% add MNW2 wells
if isfield(ssm,'MNW2') && ~isempty(ssm.MNW2)
    ssm.PNTSRC = [ssm.PNTSRC; ssm.MNW2.PNTSRC()];
end

%% add pointObj lineObj and area2Obj
if ~isfield(ssm,'PNTSRC'), ssm.PNTSRC={}; end

if isfield(ssm,'point')
    ssm.PNTSRC = [ssm.PNTSRC; ssm.line.PNTSRC(ssm.basename)];
end
if isfield(ssm,'line')
    ssm.PNTSRC = [ssm.PNTSRC; ssm.line.PNTSRC(ssm.basename)];
end
if isfield(ssm,'area')
    ssm.PNTSRC = [ssm.PNTSRC; ssm.area.PNTSRC(ssm.basename)];
end

if ~isempty(ssm.PNTSRC)

    if iscell(ssm.PNTSRC)
        ssm.PNTSRC = cell2list(ssm.PNTSRC); % SP first
    else
        ssm.PNTSRC = cell2list({ssm.PNTSRC});
    end

    if size(ssm.PNTSRC,2)<6 && ssm.NCOMP==1 || size(ssm.PNTSRC,2) < 6+ssm.NCOMP && ssm.NCOMP>1
        error(['The number of colmns is %d must be 6 or 7 if NCOMP==1 else must be 6+NCOMP=%d\n',...
                'See MT3DMS manual version 5.2, page 103 and item D8 on page 122.'],...
                 size(ssm.PNTSRC,2),6+ssm.NCOMP);
    end

    % sort, make unique and determine the first and last postion of each
    % stress period in the PNTSRC list for writing them to the SSM file
    [uniqueSP, Ilast]=  unique(ssm.PNTSRC(:,1),'last');  % unique SP in PNTSRC, with pointer to last line of each SP
    Ifirst = [0 ; Ilast(1:end-1)]+1;        % pointer to first line of each SP
    Utype  = unique(ssm.PNTSRC(:,6));       % unique ITYPES in PNTSRC

    NSS =zeros(ssm.NPER,1);        % all stress periods NSS = 0
    NSS(uniqueSP)=Ilast-Ifirst+1;  % uniqueSP have NSS~=0

    % Which option flags are on and which are off?
    optionflag (FWEL) = any (Utype == itype.WEL);  % itype well is the ITYPE number for wells (=2)
    optionflag (FDRN) = any (Utype == itype.DRN);  % itype is enumeration of class int32
    optionflag (FRIV) = any (Utype == itype.RIV);  % to see numeric itype try itype.RIV+0
    optionflag (FGHB) = any (Utype == itype.GHB);  % etc.
    optionflag (FCHD) = any (Utype == itype.CHD);
    optionflag (FMNW) = any (Utype == itype.MNW);

else
    NSS=zeros(ssm.NPER,1);
end

%%get NSS from MNW and WEL
%if isfield(ssm,'WEL');  nwel = size(vertcat(ssm.WEL.idx ),1); else nwel=0; end
%if isfield(ssm,'MNW1'); nmnw = 0;                                          end

% Setting the option flag, by writing this line the output
for iflag = 1:numel(optionflag)
    if optionflag(iflag)
        fprintf(fid,' T');
    else
        fprintf(fid,' F');
    end
end

fprintf(fid,'     FWEL FDRN FRCH FEVT FRIV FGHB FCHD FMNW FNW2 FNW3 FNW4\n');

%D2: maximum number of point sinks included in the model
fprintf(fid,'%10d     MXSS max number of sink/source points in model\n',ssm.MXSS+max(NSS));
   
%% ======================== FOR EACH PERIOD ========================

    if ssm.NCOMP==1, LAST=6; else LAST= 6 + ssm.NCOMP; end
    
    iptr = 0; % pointer tot the stress periods with non-empty PNTSRC
    
    for iPer=1:ssm.NPER

        if ssm.FRCH
            %D3 INCRCH (I10) read recharge flux for this stress period?
            if iPer>1
                for iComp=ssm.NCOMP:-1:1
                    Equal(iComp) = all(ssm.CRCH{iPer,iComp}==ssm.CRCH{iPer-1,iComp});
                end
                if all(Equal), ssm.INCRCH(iPer)=AsPrevious; end;
            end
            fprintf(fid,'%10d     INCRCH(%d)\n',ssm.INCRCH(iPer),iPer);
            if ssm.INCRCH(iPer)>=0
                for iComp=1:ssm.NCOMP
                    %D4 CRCH(NCOL,NROW) if FRCH=T and INCRCH>=0
                    warray(fid,ssm.CRCH{iPer,iComp},ssm.unit,'(10E12.3)',sprintf('CRCH{%d}',iPer));
                end
            end
        end

        if ssm.FEVT
            %D5 INCEVT (I10) read evt flux for this stress period?
            if iPer>1
                for iComp=ssm.NCOMP:-1:1
                    Equal(iComp) = all(ssm.CEVT{iPer,iComp}==ssm.CEVT{iPer-1,iComp});
                end
                if all(Equal), ssm.INCEVT(iPer)=AsPrevious; end;
            end
            fprintf(fid,'%10d     INCEVT(%d)\n',ssm.INCEVT(iPer),iPer);
            if ssm.INCEVT(iPer)>=0
                for iComp=1:ssm.NCOMP
                    %D6 CEVT(NCOL,NROW) if FEVT=T and INCEVTRCH>=0
                    warray(fid,ssm.CEVT{iPer,iComp},ssm.unit,'(10E12.3)',sprintf('CEVT{%d}',iPer));
                end
            end
        end

        %D7 NSS (I10) number of point sources whose concentrations need to
        % be specified in this stress period. That is, point sources whose
        % concentration should not be zero if water enters the model
        % through them.
        %
        % The constant-concentration condition specified in the SSM package
        % takes precedence over what is specified in the in the BTN, i.e.
        % through ICBUND. So we may override ICBUND whereever ICBUND~=-1.
        % That is point sources with ITYPE==-1 (constant head) add to what
        % already are constant heads in ICBUND (Zheng, 1999, p123).

        fprintf(fid,['%10d     NSS(%d)', ...
            '  Following lines: L R C CSS ITYPE CSSMS(1..%d) {SP %d}\n'], NSS(iPer),iPer,ssm.NCOMP,iPer);
        
        if ssm.NCOMP==1
            fmt='%10d%10d%10d%10g%10d\n';
        else
            fmt=['%10d%10d%10d%10g%10d' repmat(' %12g',[1,ssm.NCOMP]) '\n'];
        end
        if NSS(iPer)>0,
            iptr=iptr+1;
            fprintf(fid,fmt,ssm.PNTSRC(Ifirst(iptr):Ilast(iptr),2:LAST)');
        end
        
    end

    fclose(fid);
end
 
