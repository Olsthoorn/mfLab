%READ_MMF reads a MODFLOW model into Matlab (obsolete)
%
% This generic script reads files of MODFLOW, MT3D, SEAWAT, MODENSE and SWI alike
% (as far as the current implementation of the packages reaches)
% The script will be extended as more packages will be used.
% The nam read yields a cell array with extentions/package ID units and
% filenames (= the three-column nam file itself)
% Each package read yields a struct with the packag name in lower case such
% as dis, bcf, lpf, chd, riv, drn, wel, swi, cfp etc holding the package
% data.
% Note that the total of packages may become heavy in case of a large model
% Matlab matrices to be extracted from the structs are among others
%
%    DELR DELC Z KH KV SY SS PEFF IBOUND STRTHD STCONC
%
% TO 090101  090713

% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

clear variables
close all

pth=['.' filesep];

fnam=dir([pth '*.nam']);

if length(fnam)>1
    fpritnf('there are %d namfiles\n',length(fnam));
    for i=1:length(fnam)
        pfrintf('    %s\n',fnam(i).name);
    end
end

% what=3;
% switch what
%     case 0
%         testfilesmf2k=...
%         {   'bcf2ss.nam'    'ibs2k.nam'     'restest.nam'   'tc1.nam'       'tc2hufv4.nam'  'testsfr2.nam'  'twrip.nam', ...
%             'etsdrt.nam'    'l1a2k.nam'     'str.nam'       'tc1huf.nam'    'tc3.nam'       'tr2k_s3.nam', ...   
%             'example3.nam'  'l1b2k.nam'     'swtex4.nam'    'tc1obsen.nam'  'test1ss.nam'   'tvp.nam', ...
%             'fhb.nam'       'mnw1.nam'      'tc1-true.nam'  'tc2.nam'       'test1tr.nam'   'twri.nam'};
% 
%         cd('Z:\tolsthoorn On My Mac\GRWMODELS\mf2k.1_18\data');
%         namfile=testfilesmf2k{5};
%     case 1
%         cd('Z:\tolsthoorn On My Mac\GRWMODELS\MYWORK\Drainmodels\CFPM1_example');
%         cd('Z:\tolsthoorn On My Mac\GRWMODELS\MF2005\CFP_1.0\examples\CFPM1_example');
%         namfile='CFPM1.nam';
%     case 2
%         cd('Z:\tolsthoorn On My Mac\GRWMODELS\MF2005\CFP_1.0\examples\CFPM2_example');
%         namfile='CFPM2.nam';
%     case 3
%         cd('Z:\tolsthoorn On My Mac\GRWMODELS\MYWORK\AGVmodel\AGV1');
%         namfile='mf2kswi.nam';
% end


%% ===== THE NAM FILE ====================
for inam=1:length(fnam)

    nam=readNAM(fnam(inam).name,pth);

    %% ===== THE DIS-file struct ==============
    i=strmatchi('DIS',nam(:,1));
    if i, dis=readDIS(nam{i,3},pth); end

    %% ===== THE BAS FILE ====================
    i=strmatchi('BAS',nam(:,1));
    if i,
        bas.NROW=dis.NROW; bas.NCOL=dis.NCOL; bas.NLAY=dis.NLAY;
        bas=readBAS6(nam{i,3},pth,bas);
    end

    %% ===== THE BCF-file =====================
    i=strmatchi('BCF6',nam(:,1),'ErrOpt');
    if i,
        bcf.NROW=dis.NROW; bcf.NCOL=dis.NCOL; bcf.NLAY=dis.NLAY; bcf.isTran=any(dis.isTran);
        bcf.LAYCBD=dis.LAYCBD;
        bcf=readBCF6(nam{i,3},pth,bcf);
    else  % old BCF version
        i=strmatchi('BCF',nam(:,1),'ErrOpt');
        if i,
            bcf.NROW=dis.NROW; bcf.NCOL=dis.NCOL; bcf.NLAY=dis.NLAY; bcf.isTran=any(dis.isTran);
            bcf.LAYCBD=dis.LAYCBD;
            bcf=readBCF(nam{i,3},pth,bcf);
        end
    end

    %% ===== THE lpf-file ==== (use either BCF or LPF) ==============
    i=strmatchi('LPF',nam(:,1),'ErrOpt');
    if i,
        lpf.NROW=dis.NROW; lpf.NCOL=dis.NCOL; lpf.NLAY=dis.NLAY; lpf.isTran=any(dis.isTran);
        lpf.LAYCBD=dis.LAYCBD;
        lpf=readLPF(nam{i,3},pth,lpf);
    end

    %% ===== THE huf-file ===========================================
    i=strmatchi('HUF',nam(:,1),'ErrOpt');
    if i,
        huf.NROW=dis.NROW; huf.NCOL=dis.NCOL; huf.NLAY=dis.NLAY;
        huf=readHUF(nam{i,3},pth,huf);
    end
    %% ===== THE huf-file ===========================================
    i=strmatchi('HUF',nam(:,1),'ErrOpt');
    if i,
        huf.NROW=dis.NROW; huf.NCOL=dis.NCOL; huf.NLAY=dis.NLAY;
        huf=readHUF(nam{i,3},pth,huf);
    end

    %% ===== THE mult-file ===========================================
    i=strmatchi('MULT',nam(:,1),'ErrOpt');
    if i,
        mult.NROW=dis.NROW; mult.NCOL=dis.NCOL; mult.NLAY=dis.NLAY;
        mult=readMULT(nam{i,3},pth,mult);
    end


    %% ===== THE ldpa-file ===(LAYER VARIABLE DIRECTION ANISOTROPY ==
    i=strmatchi('LVDA',nam(:,1),'ErrOpt');
    if i,
        lvda.NROW=dis.NROW; lvda.NCOL=dis.NCOL; lvda.NLAY=dis.NLAY;
        lvda=readLVDA(nam{i,3},pth,lvda);
    end

    %% ===== the VDP-file for SEAWAT ================================
    i=strmatchi('VDF',nam(:,1),'ErrOpt');
    if i,
        vdf.NLAY=dis.NLAY; vdf.NROW=dis.NROW; vdf.NCOL=dis.NCOL; vdf.NPER=dis.NPER;
        vdf=readVDP(nam{i,3},pth,vdf);
    end

    %% ===== THE rch-file ============================
    i=strmatchi('RCH',nam(:,1),'ErrOpt');
    if i,
        rch.NLAY=dis.NLAY; rch.NROW=dis.NROW; rch.NCOL=dis.NCOL; rch.NPER=dis.NPER;
        rch=readRCH(nam{i,3},pth,rch);
    end

    %% s===== THE evt-file (the Evaporation package) ==============
    i=strmatchi('EVT',nam(:,1),'ErrOpt');
    if i,
        evt.NLAY=dis.NLAY; evt.NROW=dis.NROW; evt.NCOL=dis.NCOL; evt.NPER=dis.NPER;
        evt=readEVT(nam{i,3},pth,evt);
    end

    %% ===== THE WEL-file ====================================
    i=strmatchi('WEL',nam(:,1),'ErrOpt');
    if i,
        wel.NPER=dis.NPER;
    %    wel=readWEL(nam{i,3},pth,wel);
        wel=readBCN(nam{i,3},pth,wel,'WEL');

    end

    %% ===== THE GHB-file ====================================
    i=strmatchi('GHB',nam(:,1),'ErrOpt');
    if i,
        ghb.NPER=dis.NPER;
    %    ghb=readGHB(nam{i,3},pth,ghb);
        ghb=readBCN(nam{i,3},pth,ghb,'GHB');
    end

    %% ===== THE DRN-file ====================================
    i=strmatchi('DRN',nam(:,1),'ErrOpt');
    if i,
        drn.NPER=dis.NPER;
    %    drn=readDRN(nam{i,3},pth,drn);
        drn=readBCN(nam{i,3},pth,drn,'DRN');
    end

    %% ===== THE RIV-file (river package) ===============================
    i=strmatchi('RIV',nam(:,1),'ErrOpt');
    if i,
        riv.NPER=dis.NPER;
    %   riv=readRIV(nam{i,3},pth,riv);
        riv=readBCN(nam{i,3},pth,riv,'RIV');
    end

    %% ===== THE CHD-file FILE (constant head boundary package) =====
    i=strmatchi('CHD',nam(:,1),'ErrOpt');
    if i,
        chd.NPER=dis.NPER;
    %    chd=readCHD(nam{i,3},pth,chd);
        chd=readBCN(nam{i,3},pth,chd,'CHD');
    end

    %% THE BTN-file (Basic Transprot Process for MT3D)
    i=strmatchi('BTN',nam(:,1),'ErrOpt');
    if i, btn=readBTN(nam{i,3},pth); end

    %% ===== THE ADV-file (Advection process) ==========
    i=strmatchi('ADV',nam(:,1),'ErrOpt');
    if i, adv=readADV(nam{i,3},pth); end
    %% ===== THE DSP-file (Dispersion process package) =============
    i=strmatchi('DSP',nam(:,1),'ErrOpt');
    if i, dsp=readDSP(nam{i,3},pth); end

    %% ===== THE SSM-file (Source-sink mixing process package) =====
    i=strmatchi('SSM',nam(:,1),'ErrOpt');
    if i, ssm=readSSM(nam{i,3},pth); end

    %% ===== THE rct-file (chemical reaction package) ==============
    i=strmatchi('RCT',nam(:,1),'ErrOpt');
    if i, rct=readRCT(nam{i,3},pth); end

    %% ======================== SWI input ==========================
    i=strmatchi('SWI',nam(:,1),'ErrOpt');
    if i,
        swi.NROW=dis.NROW; swi.NCOL=dis.NCOL; swi.NLAY=dis.NLAY;
        swi=readSWI(nam{i,3},pth,swi);
    end

    %% ======================== CFP input ==========================
    i=strmatchi('CFP',nam(:,1),'ErrOpt');
    if i, cfp=readCFP(nam{i,3},pth); end

    %% ======================== CRCH input ==========================
    i=strmatchi('CRCH',nam(:,1),'ErrOpt');
    if i,
        crch.NPER=dis.NPER; crch.NNODES=cfp.NNODES;
        crch=readCRCH(nam{i,3},pth,crch);
    end
    %% ======================== CRCH input ==========================
    i=strmatchi('COC',nam(:,1),'ErrOpt');
    if i,
        coc=readCOC(nam{i,3},pth);
    end

end