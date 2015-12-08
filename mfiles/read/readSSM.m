function  ssm=readSSM(basename,btn)
%READSSM reads MT3DMS's source sink mixing package input file (SSM)
%
% Example:
%    ssm=readSSM(basename,btn);
%
% TO 0706030 081227 110112

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

ssm.unit=100;  % incompatibility with UD2REL (MT3DMS p97)

fid=fopen(basename,'rt');
if fid<0
    fid=fopen([basename '.SSM'],'rt');
end

fprintf('readSSM %s\n', datestr(now));

ssm.NROW=btn.NROW;
ssm.NCOL=btn.NCOL;
ssm.NLAY=btn.NLAY;
ssm.NPER=btn.NPER;
ssm.NCOMP = btn.NCOMP;
ssm.MCOMP = btn.MCOMP;
clear btn

%D0 HEADING 1+2 (<=80 chars) -- No header allowd in SSM file

%D1 FWEL FDRN FRCH FEVT FRIV FGHB FCHD FMNW (FNEW(n), n=1:4) 10I2
fprintf('reading:    FWEL FDRN FRCH FEVT FRIV FGHB FCHD FMNW (FNEW(n n=1:3)\n');
s=fgetl(fid); option_flags=sscanf(s,'%s');
ssm.FWEL=option_flags(1)=='T';
ssm.FDRN=option_flags(2)=='T';
ssm.FRCH=option_flags(3)=='T';
ssm.FEVT=option_flags(4)=='T';
ssm.FRIV=option_flags(5)=='T';
ssm.FGHB=option_flags(6)=='T';
if length(option_flags)> 6, ssm.FCHD=option_flags( 7)=='T'; end
if length(option_flags)> 7, ssm.FMNW=option_flags( 8)=='T'; end
if length(option_flags)> 8, ssm.FL09=option_flags( 9)=='T'; end
if length(option_flags)> 9, ssm.FL10=option_flags(10)=='T'; end

%% D2 maximum number of all point sinks and sources include in the flow model.
ssm.MXSS=fscanf(fid,'%d',1);
fprintf('MXSS = %d, max number of sink/source points in model\n',ssm.MXSS);
fgets(fid);
   
%% ======================== FOR EACH PERIOD ========================

for iPer=1:ssm.NPER
    
    if ssm.FRCH
        %D3 INCRCH (I10) read recharge flux for this stress period?
        fprintf('INCRCH\n');
        ssm.INCRCH(iPer) = fscanf(fid,'%d',1);
        if ssm.INCRCH(iPer)>=0
            for iComp=1:ssm.NCOMP
                %D4 CRCH(NCOL,NROW) if FRCH=T and INCRCH>=0
                ssm.CRCH{iPer,iComp}=rarray(fid,[ssm.NROW,ssm.NCOL]); %CRCH
            end
        end
    end

    if ssm.FEVT
        %D5 INCEVT (I10) read evt flux for this stress period?
        ssm.INCEVT(iPer) = fscanf(fid,'%10d',1); %ssm.INCEVT(iPer),iPer);
        if ssm.INCEVT(iPer)>=0
            for iComp=1:ssm.NCOMP
                %D6 CEVT(NCOL,NROW) if FEVT=T and INCEVTRCH>=0
                ssm.CEVT{iComp}=rarray(fid,[ssm.NROW,ssm.NCOL]); %CEVT
            end
        end
    end

    %D7 NNS (I10) max of point souces for which conc needs to be specified or
    %read for the current stress period. (MXSS is overall maximum of all sources and sinks)
    % By default, unspecified point sources are assumed zero concentration.
    % (The concentration of point sinks is always set equal
    %to the concentration of groundwater at the sink location). MT3DMS p121
    %
    % Not all sources and sinks need to be specified, only those for which
    % a concentration needs ot be specified. For constant head cells,
    % specifying the concentration means, specifying the concentration of
    % inflowing water, not that of the entire cell, which is the result of
    % mixing. This is the same as for wells and different for constand
    % concentration cells, for which the concentration of the entire cell
    % is fixed, irrespective of its size and flows in and out of it.
    %
    % So how do we specify elegantly only those point sinks for which
    % specifying the concentration of the incoming water is relevant?
    % We have to do this for each stress period.
    %
    % The specification here takes precedence over what has been specified
    % in the BTN file through ICBUND. So we may overrid ICBUND starting
    % concentrations here, if we so wish.
    ssm.NSS(iPer)=fscanf(fid,'%d',1); fgets(fid);
    
    ssm.PNTSRC=[];
    for i=1:ssm.NSS
        
        fprintf('Stress period %d, NSS= %10d point sources',iPer,ssm.NSS);
        
        % Following lines: L R C CSS ITYPE CSSMS(1..%d)\n'], NSS,iPer,ssm.NCOMP);
        
        % explore number of items on line
        p=ftell(fid); s=fgets(fid); a=scanf('%f',s); m=length(A); fseek(fid,p,'bof');
        
        pntsrc=fscanf(fid,'%f',[m,ssm.NNS])';
        ssm.PNTSRC=[ssm.PNTSRC; [iper*ones(ssm.NNS,1) pntsrc]];        
    end
end

fclose(fid);
