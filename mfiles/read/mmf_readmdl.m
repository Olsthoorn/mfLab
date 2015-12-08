%MMF_READMDL reads a MODFLOW model
%
% This generic script generates reads MODFLOW, MT3D and SEAWAT model such
% as examples coming with these programs. The result are Matlab matrices
% that defined the model. This allows easy editing of an existng model. The
% file mmf_setup is used to output a model.
%
%    DELR DELC Z KH KV SY SS PEFF IBOUND STRTHD STCONC
%
% There must also be an excel file <basename>.xls which contains al
% parameter data, stress period data and layer-specific data such as
% convertibility. See example file to understand what an how.
% TO 090101 090709


% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

%% =====RUN MMF_ADAPT ===================================================
%  MMF_ADAPT is a local m-file that may be used to the model matrices
%  created or to make a model completely
%  also in mmf_adapt you set the paths for your local model
%  TO 090101
clear variables
close all


%mmf_adapt    % local file to adapt model in current model directory

%  =====================Get started with the Seawat input ============================

clear variables
close all

%% Get the nam file to get the file names
d=dir;
% find *.nam, and use pretext as basename.

%% ===== current model ================================
fprintf('Basename current model is ''%s''\n',basename);

%XLSF=[basename '.xls'];

%% NAM Namefile INFO

% Legal packages, the list may be extented
% Specific for MODFLOW
nam.LegalMF={...
    'GLOBAL','LIST',...
    'BAS6','DIS',...
    'BCF6','LPF',...
    'RCH','EVT',...
    'WEL','RIV','GHB','DRN',...
    'UMT','CHD',...
    'PES','SEN','OBS','HOB','LMT6',...
    'PCG','OC','DATA ','DATA(BINARY)'...
    };
	
% Specific legal packages for MT3DMS
nam.LegalMT={'LIST','FTL','BTN','ADV','DSP','SSM','GCG','RCT'};

% Specific for Seawat
nam.LegalSW={nam.LegalMF{:},'VDF','BTN','ADV','DSP','SSM','GCG'};

% Specific for Salt Water Intrusion (SWI)
nam.LegalSWI={nam.LegalMF{:},'SWI'};

% STARTING WITH THE ACTUAL DATA FOR THIS MODEL =============
fprintf('Getting nam file data from %s\n',XLSF);

[Nnum,Ntxt]=xlsread(XLSF,'NAM');

% Get the data for the *.nam file from the NAM worksheet in the XSLF workbook 
nam.PCKG=Ntxt(2:end,1);     % List of packages
nam.EXT =Ntxt(2:end,3);     %   corresponding file extenstions
nam.UNIT=Nnum(:,1);         %   corresponding file units
nam.SCEN=Nnum(:,3); nam.SCEN(nam.SCEN<0)=0;  % package scenario flag
		% =0 package not inlcuded
		% >0 value is scenario number. This corresponds with data column
	    %     in the worksheets MFLOW, MT3D and SWT next to the package/parnam column
	    %     this column must be present when package is processed!!

% ====UNIT numbers for the output files (hds, ddn, bgt,zta) ==============
% These units are associate with the given extension "hds" "ddn" "bgt" "zta"
% and not with the specified package names
iunit=strmatchi('bgt',nam.EXT,1); if iunit, UBGTOUT=nam.UNIT(iunit); else UBGTOUT=0; end
iunit=strmatchi('hds',nam.EXT,1); if iunit, UHEDOUT=nam.UNIT(iunit); else UHEDOUT=0; end
iunit=strmatchi('ddn',nam.EXT,1); if iunit, UDDNOUT=nam.UNIT(iunit); else UDDNOUT=0; end
iunit=strmatchi('zta',nam.EXT,1); if iunit, UZTAOUT=nam.UNIT(iunit); else UZTAOUT=0; end
mmiunit=strmatchi('sns',nam.EXT,1); if iunit, USNSOUT=nam.UNIT(iunit); else USNSOUT=0; end

% this can be extended whenever necessary

writeNAM(basename,nam)  % GENERATE THE NAM+BAT FILES FOR MODFLOW, MT3D AND SEAWAT

%% SIMULATION INFO MODFLOW
fprintf('Getting simulation info Modflow\n');
[MFLOWparnams,MFLOWparvals]=getExcelData(XLSF,'MFLOW','Vertical');

%% SIMULATION INFO MT3D
fprintf('Getting simulation info MT3D\n');
[MT3Dparnams,MT3Dparvals]=getExcelData(XLSF,'MT3D','Vertical');

%% SIMULATION INFO SEAWAT
fprintf('Getting simulation info Seawat\n');
[SWTparnams,SWTparvals]=getExcelData(XLSF,'SWT','Vertical');

%% SIMULATION INFO SWI
fprintf('Getting simulation info SWI\n');
[SWIparnams,SWIparvals]=getExcelData(XLSF,'SWI','Vertical');

%% STRESS PERIOD INFO
fprintf('Getting stress period info\n');
[PERparnams,PERparvals]=getExcelData(XLSF,'PER','Horizontal');

%% LAYER INFO
fprintf('Getting layer info\n');
[LAYparnams,LAYparvals]=getExcelData(XLSF,'LAY','Horizontal');

% Allow fewer layer parameters if they are the same anyway
LAYCBD = LAYparvals(:,strmatchi('LAYCBD',LAYparnams));  % confining beds

NLAY=size(Z,3)-1-sum(LAYCBD(1:end-1)); %subtract the number of confining beds
[N,M]=size(LAYparvals);
if N<NLAY,  % then add layers equal to the last one specified
    LAYparvals=[LAYparvals;ones(NLAY-N,1)*LAYparvals(end,:)];
elseif N>NLAY, % then through away layers beyong NLAY
    LAYparvals(NLAY+1:N,:)=[];
end

%% ===== THE BAS FILE ====================
bas.SCEN=nam.SCEN(strmatchi('BAS6',nam.PCKG));
bas.unit=nam.UNIT(strmatchi('BAS6',nam.PCKG));
bas.ext =nam.EXT {strmatchi('BAS6',nam.PCKG)};

if bas.SCEN
    fprintf('Generating basic struct\n');
    
    bas.NLAY=NLAY; 
    bas.HNOFLO =MFLOWparvals(strmatchi('HNOFLO',MFLOWparnams),bas.SCEN);
    bas.IBOUND=IBOUND;
    bas.STRTHD=STRTHD;        % initial heads zero

    writeBAS6(basename,bas);
end


%% ===== THE DIS-file struct ==============
% STRESS PERIODS
dis.SCEN=nam.SCEN(strmatchi('DIS',nam.PCKG));
dis.unit=nam.UNIT(strmatchi('DIS',nam.PCKG));
dis.ext =nam.EXT {strmatchi('DIS',nam.PCKG)};

if dis.SCEN

    fprintf('Filling discretization struct\n');
    dis.LAYCBD = LAYparvals(:,strmatchi('LAYCBD',LAYparnams));  % confining bed below this layer
    dis.NROW=length(DELR);  % default matlab variable
    dis.NCOL=length(DELC);  % default matlab variable
    dis.NLAY=size(Z,3)-1-sum(dis.LAYCBD(1:end-1));
    dis.NPER=size(PERparvals,1);
    dis.DELR=DELR;  % column sizes
    dis.DELC=DELC;  % row sizes

    dis.LAYCBD = LAYparvals(:,strmatchi('LAYCBD',LAYparnams));  % confining bed below this layer

    % Layer elevations
    fprintf('Getting bottom elevations into\n');
    NQ3D=sum(LAYparvals(:,strmatchi('LAYCBD',LAYparnams)));

    % Specific for model of Marjolein
    dis.Z=Z;  % generated in build_model

    fprintf('Filling stress period info for discretizaton file\n'); 
    dis.PERLEN = PERparvals(:,strmatchi('PERLEN',   PERparnams));
    dis.NSTP   = PERparvals(:,strmatchi('NSTP',     PERparnams));
    dis.TSMULT = PERparvals(:,strmatchi('TSMULT',   PERparnams));
    dis.isTran = PERparvals(:,strmatchi('Transient',PERparnams));

    writeDIS(basename,dis); % write the modflow DIS file

    fprintf('Setting flag is on.\n');
end


%% ===== THE BCF-file =====================
bcf.SCEN=nam.SCEN(strmatchi('BCF6',nam.PCKG));
bcf.unit=nam.UNIT(strmatchi('BCF6',nam.PCKG));
bcf.ext =nam.EXT {strmatchi('BCF6',nam.PCKG)};

if bcf.SCEN    
    fprintf('Generating BCF struct and file\n');
    bcf.NROW=dis.NROW;
    bcf.NCOL=dis.NCOL;
    bcf.NLAY=dis.NLAY;
    bcf.DELC=dis.DELC;
    bcf.DELR=dis.DELR;
    bcf.Z=Z; % all top and bottom elevations of the model grid
    bcf.KH=KH;
    bcf.KV=KV;
    bcf.SS=SS;
    bcf.SY=SY;
    bcf.LAYCBD=dis.LAYCBD;
    %1
    %2
    bcf.IBCFCB =UBGTOUT; % flag if >0 unit nuber budget file
    bcf.HDRY   =MFLOWparvals(strmatchi('HDRY',MFLOWparnams),bcf.SCEN); % recharge option, in what layer
    bcf.IWDFLG =MFLOWparvals(strmatchi('IWDFLG',MFLOWparnams),bcf.SCEN);
    bcf.WETFCT =MFLOWparvals(strmatchi('WETFCT',MFLOWparnams),bcf.SCEN);
    bcf.IWETIT =MFLOWparvals(strmatchi('IWETIT',MFLOWparnams),bcf.SCEN);
    bcf.IHDWET =MFLOWparvals(strmatchi('IHDWET',MFLOWparnams),bcf.SCEN);

    bcf.LAYCON =LAYparvals(:,strmatchi('LAYCON' ,LAYparnams));
    bcf.LAYAVG =LAYparvals(:,strmatchi('LAYAVG' ,LAYparnams));
    bcf.TPRY   =LAYparvals(:,strmatchi('TPRY'   ,LAYparnams));
    bcf.WETDRY =LAYparvals(:,strmatchi('WETDRY' ,LAYparnams));

    bcf.isTran = PERparvals(:,strmatchi('Transient',PERparnams));

    writeBCF(basename,bcf) 
end




%% ===== THE lpf-file ==== (use either BCF or LPF) ==============
lpf.SCEN=nam.SCEN(strmatchi('LPF',nam.PCKG));
lpf.unit=nam.UNIT(strmatchi('LPF',nam.PCKG));
lpf.ext =nam.EXT {strmatchi('LPF',nam.PCKG)};

if lpf.SCEN
    fprintf('Generating LPF struct and file\n');
    lpf.NLAY=dis.NLAY;
    lpf.NROW=dis.NROW;
    lpf.NCOL=dis.NCOL;

    lpf.ILPFCB =UBGTOUT; % unit nuber budget file
    lpf.HDRY   =MFLOWparvals(strmatchi('HDRY'  ,MFLOWparnams),lpf.SCEN); % recharge option, in what layer
    lpf.NPLPF  =MFLOWparvals(strmatchi('NPLPF' ,MFLOWparnams),lpf.SCEN); % recharge option, in what layer
    
    lpf.KH=KH;
    lpf.KV=KV;
    lpf.SY=SY;
    lpf.SS=SS;

    lpf.LAYCON =LAYparvals(:,strmatchi('LAYCON' ,LAYparnams)); % Layer type 1=convertible
    lpf.LAYAVG =LAYparvals(:,strmatchi('LAYAVG' ,LAYparnams)); % Hor cond comp method (0=harmonic)
    lpf.LAYWET =LAYparvals(:,strmatchi('LAYWET' ,LAYparnams)); % Layer wettability flag
    lpf.WETDRY =LAYparvals(:,strmatchi('WETDRY' ,LAYparnams)); % Layer wetting method flag
    lpf.CHANI  =LAYparvals(:,strmatchi('CHANI'  ,LAYparnams)); % Layer wetting method flag
    lpf.LAYVKA =LAYparvals(:,strmatchi('LAYVKA' ,LAYparnams)); % 0=use kV

    lpf.isTran = PERparvals(:,strmatchi('Transient',PERparnams));

    lpf.LAYWET(lpf.LAYCON==0)=0; % LAYWET must be 0 if LAYCON is 0

    writeLPF(basename,lpf) 
end


%% ===== the VDP-file for SEAWAT ================================
vdf.SCEN=nam.SCEN(strmatchi('VDF',nam.PCKG));
vdf.unit=nam.UNIT(strmatchi('VDF',nam.PCKG));
vdf.ext =nam.EXT {strmatchi('VDF',nam.PCKG)};

if vdf.SCEN
    fprintf('Generating VDF struct and file\n');
    vdf.NLAY=dis.NLAY;
    
    % for each simulation
    %1
    vdf.MTDNCONC=SWTparvals(strmatchi('MTDNCONC',SWTparnams),vdf.SCEN); % recharge option, in what layer
    vdf.MFNADVFD=SWTparvals(strmatchi('MFNADVFD',SWTparnams),vdf.SCEN); % flag if >0 unit nuber budget file
    vdf.NSWTCPL =SWTparvals(strmatchi('NSWTCPL' ,SWTparnams),vdf.SCEN); % flag if >0 unit nuber budget file
    vdf.IWTABLE =SWTparvals(strmatchi('IWTABLE' ,SWTparnams),vdf.SCEN); % flag if >0 unit nuber budget file

    %2
    vdf.DENSEMIN =SWTparvals(strmatchi('DENSEMIN',SWTparnams),vdf.SCEN); % flag if >0 unit nuber budget file
    vdf.DENSEMAX =SWTparvals(strmatchi('DENSEMAX',SWTparnams),vdf.SCEN); % flag if >0 unit nuber budget file
 
    %3
    vdf.DNSCRIT =SWTparvals(strmatchi('DNSCRIT' ,SWTparnams),vdf.SCEN); % flag if >0 unit nuber budget file
 
    %4
    vdf.DENSEREF=SWTparvals(strmatchi('DENSEREF',SWTparnams),vdf.SCEN); % flag if >0 unit nuber budget file
    vdf.DENSESLP=SWTparvals(strmatchi('DENSESLP',SWTparnams),vdf.SCEN); % flag if >0 unit nuber budget file

    %5
    vdf.FIRSTDT =SWTparvals(strmatchi('FIRSTDT' ,SWTparnams),vdf.SCEN); % flag if >0 unit nuber budget file
    
    %% for each stress period
    %6
    vdf.INDENSE =PERparvals(:,strmatchi('INDENSE',PERparnams)); % if >0 read dense for each period
    vdf.NPER=size(vdf.INDENSE,1);
    %7
    vdf.STCONC    =STCONC;  % initial concentrations for all layers cell array   
    writeVDF(basename,vdf) 
end




%% ===== THE rch-file ============================
rch.SCEN=nam.SCEN(strmatchi('RCH',nam.PCKG));
rch.unit=nam.UNIT(strmatchi('RCH',nam.PCKG));
rch.ext =nam.EXT {strmatchi('RCH',nam.PCKG)};

if rch.SCEN
    fprintf('Generating RCH struct and file\n');
    rch.NLAY=dis.NLAY;
    rch.NROW=dis.NROW;
    rch.NCOL=dis.NCOL;
    % for each simulation
    %1
    rch.NPRCH=0;         % we will use no parameters
    %2
    rch.NRCHOP =MFLOWparvals(strmatchi('NRCHOP',MFLOWparnams),rch.SCEN); % recharge option, in what layer
    rch.IRCHCB =UBGTOUT;
    
    %3
    %4
    
    %% for each stress period
    %5
    rch.INRECH =PERparvals(:,strmatchi('INRECH',PERparnams)); % if >0 read recharge layer
    rch.INIRCH =PERparvals(:,strmatchi('INIRCH',PERparnams)); % skipped unless recharge in specific layers
    %6
    rch.RECH   =PERparvals(:,strmatchi('RECH'  ,PERparnams)); % actual rechage value (uniform)
    %7
    %8
    rch.IRCH   =PERparvals(:,strmatchi('IRCH'  ,PERparnams)); % skipped unless recharge into specific layers

    writeRCH(basename,rch) 
end
%% s===== THE evt-file (the Evaporation package) ==============
evt.SCEN=nam.SCEN(strmatchi('EVT',nam.PCKG));
evt.unit=nam.UNIT(strmatchi('EVT',nam.PCKG));
evt.ext =nam.EXT {strmatchi('EVT',nam.PCKG)};

if evt.SCEN
    fprintf('Generating EVT struct and file\n');
    evt.NLAY=dis.NLAY;
    evt.NROW=dis.NROW;
    evt.NCOL=dis.NCOL;
    %1
    %2
    evt.NEVTOP =MFLOWparvals(strmatchi('NEVTOP',MFLOWparnams),evt.SCEN);
    evt.IEVTCB =UBGTOUT;
    %3
    %4
    %5
    evt.INSURF =PERparvals(:,strmatchi('INSURF',PERparnams));
    evt.INEVTR =PERparvals(:,strmatchi('INEVTR',PERparnams));
    evt.INEXDP =PERparvals(:,strmatchi('INEXDP',PERparnams));
    evt.INIEVT =PERparvals(:,strmatchi('INIEVT',PERparnams));
    %6
    evt.SURF =cell(dis.NPER,1);
    for iL=1:dis.NPER
    varname=PERparvals{iL,strmatchi('SURF',PERparnams)};
        %eval(['load ',varname]);
        eval(['evt.SURF{iL} =',varname,';']);
    end
    %7
    evt.EVTR =PERparvals(:,strmatchi('EVTR',PERparnams));
    evt.evap =evap;
    %8
    %9
    evt.EXDP =PERparvals(:,strmatchi('EXDP',PERparnams));
    %10
    evt.IEVT =PERparvals(:,strmatchi('IEVT',PERparnams));

    writeEVT(basename,evt,evtp) 
end


%% ===== THE WEL-file ====================================
wel.SCEN=nam.SCEN(strmatchi('WEL',nam.PCKG));
wel.unit=nam.UNIT(strmatchi('WEL',nam.PCKG));
wel.ext =nam.EXT {strmatchi('WEL',nam.PCKG)};

if wel.SCEN
    wel.IWELCB=UBGTOUT;  % budget file unit
    
    fprintf('Getting well info\n');
    [wel.parnams,wel.parvals]=getExcelData(XLSF,'WEL','Horizontal');
    
    PERCOL=strmatchi('PERIOD',wel.parnams);
    
    % ITMP is number of wells in each period, 0=reuse wells of last period
    % the number of periods is already known from reading per sheet
    wel.ITMP=zeros(dis.NPER,1); % set all to reuse, no wells until the first specified period is encoutnered
    
    wel.parvals=sortrows(wel.parvals);  % make sure periods are ascending

    for iP=1:dis.NPER
        wel.ITMP(iP)=sum(wel.parvals(:,PERCOL)==iP);  % count wells in this period
    end
    
    writeWEL(basename,wel)
end



%% ===== THE DRN-file ====================================
drn.SCEN=nam.SCEN(strmatchi('DRN',nam.PCKG));
drn.unit=nam.UNIT(strmatchi('DRN',nam.PCKG));
drn.ext =nam.EXT {strmatchi('DRN',nam.PCKG)};

if drn.SCEN
    fprintf('Generating DRN struct and file\n');
    [drn.values,drn.headers]=xlsread(XLSF,'DRN');
    drn.IDRNCB =UBGTOUT;
    drn.ITMP   =PERparvals(:,strmatchi('ITMPD'  ,PERparnams));
    % the ones(size(drains,1)) adds period number 1
    drn.values=[drn.values;[ones(size(drains,1),1), drains]];

    writeDRN(basename,drn,drnp)
end


%% ===== THE CHD-file FILE (constant head boundary package) =====
chd.SCEN=nam.SCEN(strmatchi('CHD',nam.PCKG));
%chd.unit=nam.UNIT(strmatchi('CHD',nam.PCKG)); % No separate CHD data to budget file
chd.ext =nam.EXT {strmatchi('CHD',nam.PCKG)};

if chd.SCEN
    fprintf('Getting CHD info\n');
    [chd.parnams,chd.parvals]=getExcelData(XLSF,'CHD','Horizontal');
    
    PERCOL=strmatchi('PERIOD',chd.parnams);
    
    % ITMP is number of wells in each period, 0=reuse wells of last period
    % the number of periods is already known from reading per sheet
    chd.ITMP=zeros(dis.NPER,1); % set all to reuse, no wells until the first specified period is encoutnered
    
    chd.parvals=sortrows(chd.parvals);  % make sure periods are ascending

    for iP=1:dis.NPER
        chd.ITMP(iP)=sum(chd.parvals(:,PERCOL)==iP);  % count wells in this period
    end
    
    writeCHD(basename,chd)
end


%% ===== THE RIV-file (river package) ===============================
riv.SCEN=nam.SCEN(strmatchi('RIV',nam.PCKG));
riv.unit=nam.UNIT(strmatchi('RIV',nam.PCKG));
riv.ext =nam.EXT {strmatchi('RIV',nam.PCKG)};

if riv.SCEN
    fprintf('Generating RIV struct and file\n');
    [riv.values,riv.headers]=xlsread(XLSF,'RIV');
    riv.IRIVCB =UBGTOUT;
    riv.ITMP   =PERparvals(:,strmatchi('ITMPD'  ,PERparnams),riv.SCEN);
    riv.IMTP   =PERparvals(:,strmatchi('ITMPR'  ,PERparnams),riv.SCEN);
    % the ones(size(drains,1)) adds period number 1
    riv.values=[riv.values;[ones(size(rivresm,1),1), rivresm]];
    writeRIV(basename,riv)
end


%% ===== THE GHB-file (general head boundaries package) ===============================
ghb.SCEN=nam.SCEN(strmatchi('GHB',nam.PCKG));
ghb.unit=nam.UNIT(strmatchi('GHB',nam.PCKG));
ghb.ext =nam.EXT {strmatchi('GHB',nam.PCKG)};

if ghb.SCEN
    fprintf('Generating GHB struct and file\n');
    [ghb.values,ghb.headers]=xlsread(XLSF,'GHB');
    ghb.IGHBCB =UBGTOUT;
    ghb.NPER=dis.NPER;
    % the ones(size(drains,1)) adds period number 1
    %ghb.values=[ghb.values;[ones(size(ghbresm,1),1), ghbresm]];
    writeGHB(basename,ghb)
end



%% ==== THE OC-file (output control package) =======================
oc.SCEN=nam.SCEN(strmatchi('OC',nam.PCKG));
oc.unit=nam.UNIT(strmatchi('OC',nam.PCKG));
oc.ext =nam.EXT {strmatchi('OC',nam.PCKG)};

if oc.SCEN
    % for each simulation
    fprintf('Generating OC struct and file\n');
    oc.IHEDFM=MFLOWparvals(strmatchi('IHEDFM',MFLOWparnams),oc.SCEN);
    oc.IDDNFM=MFLOWparvals(strmatchi('IDDNFM',MFLOWparnams),oc.SCEN);
    oc.IHEDUN=UHEDOUT;  % read from namfile not from spreadsheet
    oc.IDDNUN=UDDNOUT;  % read from namfile not from spreadsheet

    % for each stres period
    oc.INCODE =PERparvals(:,strmatchi('INCODE',PERparnams));
    oc.IHDDFL =PERparvals(:,strmatchi('IHDDFL',PERparnams));
    oc.IBUDFL =PERparvals(:,strmatchi('IBUDFL',PERparnams));
    oc.ICBCFL =PERparvals(:,strmatchi('ICBCFL',PERparnams));
    oc.NPER=size(oc.INCODE,1);
    oc.NSTP=dis.NSTP;
    
    oc.Hdpr =PERparvals(:,strmatchi('Hdpr',PERparnams));
    oc.Ddpr =PERparvals(:,strmatchi('Ddpr',PERparnams));
    oc.Hdsv =PERparvals(:,strmatchi('Hdsv',PERparnams));
    oc.Ddsv =PERparvals(:,strmatchi('Ddsv',PERparnams));

    writeOC(basename,oc);
end




%% ===== THE pcg-file (P Conjugate solver Package) =====
pcg.SCEN=nam.SCEN(strmatchi('PCG',nam.PCKG));
pcg.unit=nam.UNIT(strmatchi('PCG',nam.PCKG));
pcg.ext =nam.EXT {strmatchi('PCG',nam.PCKG)};

if pcg.SCEN
    fprintf('Generating PCG struct and file\n');
    pcg.MXITER=MFLOWparvals(strmatchi('MXITER',MFLOWparnams),pcg.SCEN);
    pcg.ITER1 =MFLOWparvals(strmatchi('ITER1' ,MFLOWparnams),pcg.SCEN);
    pcg.NPCOND=MFLOWparvals(strmatchi('NPCOND',MFLOWparnams),pcg.SCEN);
    pcg.HCLOSE=MFLOWparvals(strmatchi('HCLOSE',MFLOWparnams),pcg.SCEN);
    pcg.RCLOSE=MFLOWparvals(strmatchi('RCLOSE',MFLOWparnams),pcg.SCEN);
    pcg.RELAX =MFLOWparvals(strmatchi('RELAX' ,MFLOWparnams),pcg.SCEN);
    pcg.NBPOL =MFLOWparvals(strmatchi('NBPOL' ,MFLOWparnams),pcg.SCEN);
    pcg.IPRPCG=MFLOWparvals(strmatchi('IPRPCG',MFLOWparnams),pcg.SCEN);
    pcg.MUTPCG=MFLOWparvals(strmatchi('MUTPCG',MFLOWparnams),pcg.SCEN);
    pcg.DAMP  =MFLOWparvals(strmatchi('DAMP'  ,MFLOWparnams),pcg.SCEN);

    writePCG(basename,pcg);
end

%% ===== THE LMT6 file (Link file generation LMT6) =====
lmt.SCEN=nam.SCEN(strmatchi('LMT',nam.PCKG));
lmt.unit=nam.UNIT(strmatchi('LMT',nam.PCKG));  % must not be used !!
lmt.ext =nam.EXT {strmatchi('LMT',nam.PCKG)};

if lmt.SCEN
    lmt.OUTPUT_FILE_NAME='';
    lmt.OUTPUT_FILE_HEADER='STANDARD';
    lmt.OUTPUT_FILE_FORMAT='UNFORMATTED';

    writeLMT(basename,lmt);
end


%% ===== END OF INPUT FOR THE FLOW PROCESs ==================

%% ==== MT3DSM input ==================================
% partly from the MT3D sheet in the spreadsheet


%% THE BTN-file (Basic Transprot Process for MT3D)
btn.SCEN=nam.SCEN(strmatchi('BTN',nam.PCKG));
btn.unit=nam.UNIT(strmatchi('BTN',nam.PCKG));
btn.ext =nam.EXT {strmatchi('BTN',nam.PCKG)};

if btn.SCEN    
    NCOMP=1;  % # total of chemcial species involved
    MCOMP=1;  % # total of mobile   chemical species
    
    fprintf('Generating Basic Transport Process struct\n');
    %3
    btn.NLAY=dis.NLAY;   % # layers
    btn.NROW=dis.NROW;   % # rows
    btn.NCOL=dis.NCOL;   % # columns
    btn.NPER=size(PERparvals,1);   % # stress periods
    btn.NCOMP=NCOMP; % # chemical species
    btn.MCOMP=MCOMP; % # mobile mobile species
    %4
    btn.TUNIT='DAY';
    btn.LUNIT='M';
    btn.MUNIT='KG';
    
    %5
    %Packge-use flags: (ADV DSP SSM RCT GCG XXX XXX XXX XXX XXX)';
    adv.SCEN=nam.SCEN(strmatchi('ADV',nam.PCKG));  % if true then active
    dsp.SCEN=nam.SCEN(strmatchi('DSP',nam.PCKG));  % etc
    ssm.SCEN=nam.SCEN(strmatchi('SSM',nam.PCKG));
    rct.SCEN=nam.SCEN(strmatchi('RCT',nam.PCKG));
    gcg.SCEN=nam.SCEN(strmatchi('GCG',nam.PCKG));
    % active flags then become
    btn.TRNOP=[adv.SCEN,dsp.SCEN,ssm.SCEN,rct.SCEN,gcg.SCEN,0,0,0,0,0];
    
    %6
    btn.LAYCON=LAYparvals(:,strmatchi('LAYCON',LAYparnams));
    %7
    btn.DELR=dis.DELR;
    %8
    btn.DELC=dis.DELC;
    %9 %10
    btn.Z=Z;
    %11
    btn.PRSITY=PEFF;
    %12
    btn.ICBUND=ICBUND;  % 0 inactive, 1 active, -1 fixed concentration 
    %13
    btn.STCONC=STCONC;
    %14
    btn.CINACT=MT3Dparvals(strmatchi('CINACT',MT3Dparnams),btn.SCEN);
    btn.THKMIN=MT3Dparvals(strmatchi('THKMIN',MT3Dparnams),btn.SCEN);
    %15
    btn.IFMTCN=MT3Dparvals(strmatchi('IFMTCN',MT3Dparnams),btn.SCEN);
    btn.IFMTNP=MT3Dparvals(strmatchi('IFMTNP',MT3Dparnams),btn.SCEN);
    btn.IFMTRF=MT3Dparvals(strmatchi('IFMTRF',MT3Dparnams),btn.SCEN);
    btn.IFMTDP=MT3Dparvals(strmatchi('IFMTDP',MT3Dparnams),btn.SCEN);
    btn.SAVUCN=MT3Dparvals(strmatchi('SAVUCN',MT3Dparnams),btn.SCEN);
    %16 %17
    btn.NPRS  =MT3Dparvals(strmatchi('NPRS',MT3Dparnams),btn.SCEN);
    btn.TIMPRS=MT3Dparvals(strmatchi('TIMPRS',MT3Dparnams),:); % record with times to print output    
    %18
    btn.NOBS  =MT3Dparvals(strmatchi('NOBS',  MT3Dparnams),btn.SCEN);
    btn.NPROBS=MT3Dparvals(strmatchi('NPROBS',MT3Dparnams),btn.SCEN);
     %19  skipped, are the locations of the observation points
    %% OBSERVATION POINTS FOR CONCENTRATION
    if btn.NOBS>0
        fprintf('Getting observatoin points\n');
        [btn.OBSparnams,btn.OBS_LRC]=getExcelData(XLSF,'BTNOBS','Horizontal');
    end
    %20
    btn.CHKMAS=MT3Dparvals(strmatchi('CHKMAS',MT3Dparnams),btn.SCEN);
    btn.NPRMAS=MT3Dparvals(strmatchi('NPRMAS',MT3Dparnams),btn.SCEN);
    %21
    btn.PERLEN = dis.PERLEN;
    btn.NSTP   = dis.NSTP;
    btn.TSMULT = dis.TSMULT;
    %22 if TSMULT<0 skip
    %23
    btn.DT0    =PERparvals(:,strmatchi('DT0'    ,PERparnams)); % i.e. let the model choose
    btn.MXSTRN =PERparvals(:,strmatchi('MXSTRN' ,PERparnams)); % maximum number of transport steps within one flow time step
    btn.TTSMULT=PERparvals(:,strmatchi('TTSMULT',PERparnams)); % succuessive transport time step multiplier
    btn.TTSMAX =PERparvals(:,strmatchi('TTSMAX' ,PERparnams)); %  use default, set no maximum
    
    writeBTN(basename,btn);
end




%% ===== THE ADV-file (Advection process) ==========
adv.SCEN=nam.SCEN(strmatchi('ADV',nam.PCKG));
adv.unit=nam.UNIT(strmatchi('ADV',nam.PCKG));
adv.ext =nam.EXT {strmatchi('ADV',nam.PCKG)};

if adv.SCEN   
    fprintf('Generating Advection Process struct\n');
    %3
    adv.MIXELM=MT3Dparvals(strmatchi('MIXELM' ,MT3Dparnams),adv.SCEN);   % standard finite difference method
    adv.PERCEL=MT3Dparvals(strmatchi('PERCEL' ,MT3Dparnams),adv.SCEN);
    adv.MXPART=MT3Dparvals(strmatchi('MXPART' ,MT3Dparnams),adv.SCEN);
    adv.NADVFD=MT3Dparvals(strmatchi('NADVFD' ,MT3Dparnams),adv.SCEN);
    %2
    adv.ITRACK=MT3Dparvals(strmatchi('ITRACK' ,MT3Dparnams),adv.SCEN);   % particle tracking algorith flag 3=good compromise
    adv.WD    =MT3Dparvals(strmatchi('WD'     ,MT3Dparnams),adv.SCEN);   % concentration weighting factor, 05 should be adequat
    %3
    adv.DCEPS =MT3Dparvals(strmatchi('DCEPS'  ,MT3Dparnams),adv.SCEN);
    adv.NPLANE=MT3Dparvals(strmatchi('NPLANE' ,MT3Dparnams),adv.SCEN);
    adv.NPL   =MT3Dparvals(strmatchi('NPL'    ,MT3Dparnams),adv.SCEN);
    adv.NPH   =MT3Dparvals(strmatchi('NPH'    ,MT3Dparnams),adv.SCEN);
    adv.NPMIN =MT3Dparvals(strmatchi('NPMIN'  ,MT3Dparnams),adv.SCEN);
    adv.NPMAX =MT3Dparvals(strmatchi('NPMAX'  ,MT3Dparnams),adv.SCEN);
    %4
    adv.INTERP =MT3Dparvals(strmatchi('INTERP',MT3Dparnams),adv.SCEN);
    adv.NLSINK =MT3Dparvals(strmatchi('NLSINK',MT3Dparnams),adv.SCEN);
    adv.NPSINK =MT3Dparvals(strmatchi('NPSINK',MT3Dparnams),adv.SCEN);
    
    adv.DCHMOC =MT3Dparvals(strmatchi('DCHMOC',MT3Dparnams),adv.SCEN);
    writeADV(basename,adv);
end




%% ===== THE DSP-file (Dispersion process package) =============
dsp.SCEN=nam.SCEN(strmatchi('DSP',nam.PCKG));
dsp.unit=nam.UNIT(strmatchi('DSP',nam.PCKG));
dsp.ext =nam.EXT {strmatchi('DSP',nam.PCKG)};

if dsp.SCEN
   fprintf('Generating Dispersion Process struct\n');
    %3
    dsp.NLAY=dis.NLAY;
    dsp.NROW=dis.NROW;
    dsp.NCOL=dis.NCOL;
    
    dsp.AL    =LAYparvals(:,strmatchi('AL'    ,LAYparnams));
    dsp.TRPT  =LAYparvals(:,strmatchi('TRPT'  ,LAYparnams));  % aTH/aL
    dsp.TRPV  =LAYparvals(:,strmatchi('TRPV'  ,LAYparnams));  % aTV/aL
    dsp.DMCOEF=LAYparvals(:,strmatchi('DMCOEF',LAYparnams));  % effctief molecular diffusion coefficient
    writeDSP(basename,dsp);
end




%% ===== THE SSM-file (Source-sink mixing process package) =====
ssm.SCEN=nam.SCEN(strmatchi('SSM',nam.PCKG));
ssm.unit=nam.UNIT(strmatchi('SSM',nam.PCKG));
ssm.ext =nam.EXT {strmatchi('SSM',nam.PCKG)};

if ssm.SCEN
    fprintf('Generating Source-Sink Mixing Process struct\n');
    %D1 FWEL FDRN FRCH FEVT FRIV FGHB (FNEW(n), n=1:4) 10I2
    % option flags (T F) for WEL, DRN etc packages
    % FNEW (n=1:4) future option flags
    %3 FLAGS for packages used by FLOW process
    %conversion from  0 1 to 'F' 'T' is done in writeSSM
    ssm.FWEL  =wel.SCEN;
    ssm.FDRN  =drn.SCEN;
    ssm.FRCH  =rch.SCEN;
    ssm.FEVT  =evt.SCEN;
    ssm.FRIV  =riv.SCEN;
    ssm.FGHB  =ghb.SCEN;
    ssm.FNW1  =0;
    ssm.FNW2  =0;
    ssm.FNW3  =0;
    ssm.FNW4  =0;
    
    ssm.INCRCH =PERparvals(:,strmatchi('INCRCH',PERparnams));
    ssm.CRCH   =PERparvals(:,strmatchi('CRCH_1',  PERparnams));
    ssm.INCEVT =PERparvals(:,strmatchi('INCEVT',PERparnams));
    ssm.CEVT   =PERparvals(:,strmatchi('CEVT_1',  PERparnams));
    ssm.NPER=size(PERparvals,1);
    ssm.NROW=dis.NROW;
    ssm.NCOL=dis.NCOL;
    ssm.NLAY=dis.NLAY;
    ssm.NCOMP=btn.NCOMP;
    ssm.IBOUND=IBOUND;
    ssm.ICBUND=ICBUND;
    ssm.STCONC=STCONC;
    if ssm.FWEL, ssm.well=wel.parvals; end

    try
        fprintf('Getting point sourses of different types, see MT3MS manual, D8, p122\n');
        [PNTSRC.parnams,PNTSRC.parvals]=getExcelData(XLSF,'PNTSRC','Horizontal');
        ssm.pntsrc=PNTSRC.parvals;
    catch
        fprintf('No PNTSRC worksheet in xls file %s',XLSF);
        ssm.pntsrc=[];
    end

    writeSSM(basename,ssm);
end


%% ===== THE GCG-file (Generalize conjugate gradient solver package)
gcg.SCEN=nam.SCEN(strmatchi('GCG',nam.PCKG));
gcg.unit=nam.UNIT(strmatchi('GCG',nam.PCKG));
gcg.ext =nam.EXT {strmatchi('GCG',nam.PCKG)};

if gcg.SCEN
    fprintf('Generating Generalized Conjugate Gradient Solver Process struct\n');
    %3
    gcg.MXITER=MT3Dparvals(strmatchi('MXITER',MT3Dparnams),gcg.SCEN);
    gcg.ITER1 =MT3Dparvals(strmatchi('ITER1' ,MT3Dparnams),gcg.SCEN);
    gcg.ISOLVE=MT3Dparvals(strmatchi('ISOLVE',MT3Dparnams),gcg.SCEN);
    gcg.NCRS  =MT3Dparvals(strmatchi('NCRS'  ,MT3Dparnams),gcg.SCEN);
    gcg.ACCL  =MT3Dparvals(strmatchi('ACCL'  ,MT3Dparnams),gcg.SCEN);
    gcg.CCLOSE=MT3Dparvals(strmatchi('CCLOSE',MT3Dparnams),gcg.SCEN);
    gcg.IPRGCG=MT3Dparvals(strmatchi('IPRGCG',MT3Dparnams),gcg.SCEN);

    writeGCG(basename,gcg);
end



%% ===== THE rct-file (chemical reaction package) =======
rct.SCEN=nam.SCEN(strmatchi('RCT',nam.PCKG));
rct.unit=nam.UNIT(strmatchi('RCT',nam.PCKG));
rct.ext =nam.EXT {strmatchi('RCT',nam.PCKG)};

if rct.SCEN
    rct.NCOL=dis.NCOL;
    rct.NROW=dis.NROW;
    rct.NLAY=dis.NLAY;
        
    fprintf('Generating Chemical Reaction Process struct\n');
    %E1 ISOTHM IREACT IRCTOP IGETSC
    %  ISOTHM sorption type
    %    0=no sorption
    %    1=linear sorption
    %    2=freundlich isotherm
    %    3=langmuir
    %    4=first order kinetic sorption (non-equilibrium)
    %    5=dual domain mass transfer (without sorption)
    %    6=dual domain mass transfer (with sorption)
    rct.ISOTHM=MT3Dparvals(strmatchi('ISOTHM',MT3Dparnams),rct.SCEN);

    % IREACT type of kinetic rate reaction flag
    %    0=no kinetic rate reaction
    %    1=first-order kinetic rate reaction
    %      chemical reaction simulations need add-on package
    rct.IREACT=MT3Dparvals(strmatchi('IREACT',MT3Dparnams),rct.SCEN);
  
    % IRCTOP reaction variable entry method flag
    %  >=2 all variables are entered as 3D array (using RARRAY)
    %  < 2 all variables are entered as 1D with one value per layer
    rct.IRCTOP=MT3Dparvals(strmatchi('IRCTOP',MT3Dparnams),rct.SCEN);
    
    % IGETSC initial conc for the adsorbed phase reading flag (ISOTHM=4,5,6)
     rct.IGETSC=MT3Dparvals(strmatchi('IGETSC',MT3Dparnams),rct.SCEN);
     
     rct.NLAY=size(LAYparvals,1);
     
     rct.RHOB   =LAYparvals(:,strmatchi('RHOB',   LAYparnams));
     rct.PRSITY2=LAYparvals(:,strmatchi('PRSITY2',LAYparnams));
    
     rct.NCOMP=btn.NCOMP;
     
     rct.SRCONC=zeros(rct.NLAY,rct.NCOMP);
     rct.SP1   =zeros(rct.NLAY,rct.NCOMP);
     rct.SP2   =zeros(rct.NLAY,rct.NCOMP);
     rct.RC1   =zeros(rct.NLAY,rct.NCOMP);
     rct.RC2   =zeros(rct.NLAY,rct.NCOMP);
     for iCOMP=1:rct.NCOMP
         s=sprintf('_%d',iCOMP);
         rct.SRCONC(:,iCOMP) =LAYparvals(:,strmatchi(['SRCONC',s], LAYparnams));
         rct.SP1(:,iCOMP)    =LAYparvals(:,strmatchi(['SP1',   s], LAYparnams));
         rct.SP2(:,iCOMP)    =LAYparvals(:,strmatchi(['SP2',   s], LAYparnams));
         rct.RC1(:,iCOMP)    =LAYparvals(:,strmatchi(['RC1',   s], LAYparnams));
         rct.RC2(:,iCOMP)    =LAYparvals(:,strmatchi(['RC2',   s], LAYparnams));
     end
   
    writeRCT(basename,rct);
end

%% ===== END of MT3DSM input ==============================
%% ======================== SWI input ============================
swi.SCEN=nam.SCEN(strmatchi('SWI',nam.PCKG));
swi.unit=nam.UNIT(strmatchi('SWI',nam.PCKG));
swi.ext =nam.EXT {strmatchi('SWI',nam.PCKG)};

if swi.SCEN

    fprintf('Generating Salt Water Intrusion struct\n');
    [INTERPHparnams INTERPHparvals]=getExcelData(XLSF,'INTERPH','Horizontal');
    [ZONEparnams ZONEparvals]=getExcelData(XLSF,'ZONE','Horizontal');
    
    swi.ZONES = size(ZONEparvals,1);
    swi.Z=Z;
    swi.LAYCBD=dis.LAYCBD;
    swi.NPLN = swi.ZONES-1;
    swi.ISTRAT   = SWIparvals(strmatchi('ISTRAT',     SWIparnams),:);
    swi.ISWIZT = SWIparvals(strmatchi('ISWIZT',   SWIparnams),:);
    swi.NPRN = SWIparvals(strmatchi('NPRN',SWIparnams),:);
    swi.TOESLOPE   = SWIparvals(strmatchi('TOESLOPE',     SWIparnams),:);
    swi.TIPSLOPE = SWIparvals(strmatchi('TIPSLOPE',   SWIparnams),:);
    swi.ZETAMIN = SWIparvals(strmatchi('ZETAMIN',SWIparnams),:);
    swi.DELZETA = SWIparvals(strmatchi('DELZETA',SWIparnams),:);
    
    if swi.ISTRAT
        swi.NU=ZONEparvals(:,strmatchi('NU', ZONEparnams));
    else
        swi.NU=INTERPHparvals(:,strmatchi('NU', INTERPHparnams));
    end
    
    swi.IFace=IFace;
    swi.SSZ=LAYparvals(:,strmatchi('SSZ',LAYparnams));
    swi.ISOURCE=ISOURCE;
    swi.NLAY=NLAY;
    
    writeSWI(basename,swi);
end


%% ============ End of SWI input ============

fprintf('Ready ... %s!! \n\n\n',datestr(now));


%% ========= RUN Seawat MT3DMS or mf2k =====================

%dos 'swt_v4.bat'
%dos 'mt3dms5b.bat'
%dos 'mf2k.bat'

%dos('copy D:\dharomonteagudo\Desktop\GRWMODELS\swiex\ex1\swiex1.swi Boxswi.swi')



%dos 'mf2kswi.bat'

%% ==== use mmf_analyze to analyze the output

%mmf_analyze

% mmf_analyze is a working m-file to be adapted as to
% optimize the performance or vizualization
% you can use a copy of mmfm_analyze at every working directory
% matlab will used the one on the current directoty.
% you can always check which script is used by typing
% which mmf_analyze
