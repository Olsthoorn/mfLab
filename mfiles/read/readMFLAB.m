%READMFLAB reads a MODFLOW model into mfLab
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

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fprintf('Reading model input data readMFLAB\n');

pth=['./Hanik' filesep];

%fnam=dir([pth 'HanikMT3D.nam']);
fnam=dir([pth 'Hanik.nam']);

if ~isempty(fnam)
    fprintf('there are %d namfiles\n',length(fnam));
    for i=1:length(fnam)
        fprintf('    %s\n',fnam(i).name);
    end
else
    error('No .nam files found, don''t know what to do and quit!\n');
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

    namfile=fnam(inam).name;     nam=readNAM(namfile,pth);
    
        %% ===== THE DIS-file struct ==============
    try
        i=strmatchi('DIS',nam(:,1));
        dis=readDIS([pth nam{i,3}]);
    catch ERR
        fprintf('No DIS in %s\n',namfile);
        fprintf('identifyer: %s\n',ERR.identifier);
        fprintf('message:    %s\n',ERR.message);
    end

        %% ===== THE BCF-file =====================
    try
        i=strmatchi('BCF',nam(:,1)); % old BCF version
        bcf.NROW   = dis.NROW;
        bcf.NCOL   = dis.NCOL;
        bcf.NLAY   = dis.NLAY;
        bcf.LAYCBD = dis.LAYCBD;
        bdf.NPER   = dis.NPER;
        bcf=readBCF([pth nam{i,3}],bcf);
    catch ERR
        fprintf('No BCF in %s\n',namfile);
        fprintf('identifyer: %s\n',ERR.identifier);
        fprintf('message:    %s\n',ERR.message);
    end
    
    %% ===== THE BCF-file =====================
    try
        i=strmatchi('BCF6',nam(:,1),'exact');
        bcf.NROW   = dis.NROW;
        bcf.NCOL   = dis.NCOL;
        bcf.NLAY   = dis.NLAY;
        bcf.LAYCBD = dis.LAYCBD;
        bcf.NPER   = dis.NPER;
        bcf=readBCF6([pth nam{i,3}],bcf);
    catch ERR
        fprintf('No BCF6 in %s\n',namfile);
        fprintf('identifyer: %s\n',ERR.identifier);
        fprintf('message:    %s\n',ERR.message);
    end
    
    %% ===== THE BAS6 FILE ====================
    try
        i=strmatchi('BAS6',nam(:,1),'exact');
        bas.NROW=dis.NROW;
        bas.NCOL=dis.NCOL;
        bas.NLAY=dis.NLAY;
        bas=readBAS6([pth nam{i,3}],bas);
    catch ERR
        fprintf('No BAS6 in %s\n',namfile);
        fprintf('identifyer: %s\n',ERR.identifier);
        fprintf('message:    %s\n',ERR.message);
    end
    
    try
        i=strmatchi('BAS',nam(:,1),'exact');
        bas=readBAS([pth nam{i,3}]);
    catch ERR
        fprintf('No BAS in %s\n',namfile);
        fprintf('identifyer: %s\n',ERR.identifier);
        fprintf('message:    %s\n',ERR.message);
    end
    
    %% ===== THE lpf-file ==== (use either BCF or LPF) ==============
    try
        i=strmatchi('LPF',nam(:,1));
        lpf.NROW=dis.NROW; lpf.NCOL=dis.NCOL; lpf.NLAY=dis.NLAY; lpf.isTran=any(dis.isTran);
        lpf.LAYCBD=dis.LAYCBD;
        lpf=readLPF([pth nam{i,3}],lpf);
    catch ERR
        fprintf('No LPF in %s\n',namfile);
        fprintf('identifyer: %s\n',ERR.identifier);
        fprintf('message:    %s\n',ERR.message);
    end

    %% ===== THE huf-file ===========================================
    try
        i=strmatchi('HUF',nam(:,1));
        huf.NROW=dis.NROW; huf.NCOL=dis.NCOL; huf.NLAY=dis.NLAY;
        huf=readHUF([pth nam{i,3}],huf);
    catch ERR
        fprintf('No HUF in %s\n',namfile);
        fprintf('identifyer: %s\n',ERR.identifier);
        fprintf('message:    %s\n',ERR.message);
    end

    %% ===== THE mult-file ===========================================
    try
        i=strmatchi('MULT',nam(:,1));
        mult.NROW=dis.NROW; mult.NCOL=dis.NCOL; mult.NLAY=dis.NLAY;
        mult=readMULT([pth nam{i,3}],mult);
    catch ERR
        fprintf('No MULT in %s\n',namfile);
        fprintf('identifyer: %s\n',ERR.identifier);
        fprintf('message:    %s\n',ERR.message);
    end


    %% ===== THE ldpa-file ===(LAYER VARIABLE DIRECTION ANISOTROPY ==
    try
        i=strmatchi('LVDA',nam(:,1));
        lvda.NROW=dis.NROW; lvda.NCOL=dis.NCOL; lvda.NLAY=dis.NLAY;
        lvda=readLVDA([pth nam{i,3}],lvda);
    catch ERR
        fprintf('No LVDA in %s\n',namfile);
        fprintf('identifyer: %s\n',ERR.identifier);
        fprintf('message:    %s\n',ERR.message);
   end

    %% ===== the VDF-file for SEAWAT ================================
    try
        i=strmatchi('VDF',nam(:,1));
        vdf.NLAY=dis.NLAY; vdf.NROW=dis.NROW; vdf.NCOL=dis.NCOL; vdf.NPER=dis.NPER;
        vdf=readVDP([pth nam{i,3}],vdf);
    catch ERR
        fprintf('No VDF in %s\n',namfile);
        fprintf('identifyer: %s\n',ERR.identifier);
        fprintf('message:    %s\n',ERR.message);
     end

    %% ===== THE rch-file ============================
    try
        i=strmatchi('RCH',nam(:,1));
        rch.NLAY=dis.NLAY; rch.NROW=dis.NROW; rch.NCOL=dis.NCOL; rch.NPER=dis.NPER;
        rch=readRCH([pth nam{i,3}],rch);
    catch ERR
        fprintf('No RCH in %s\n',namfile);
        fprintf('identifyer: %s\n',ERR.identifier);
        fprintf('message:    %s\n',ERR.message);
     end

    %% s===== THE evt-file (the Evaporation package) ==============
    try
        i=strmatchi('EVT',nam(:,1));
        evt.NLAY=dis.NLAY; evt.NROW=dis.NROW; evt.NCOL=dis.NCOL; evt.NPER=dis.NPER;
        evt=readEVT([pth nam{i,3}],evt);
    catch ERR
        fprintf('No EVT in %s\n',namfile);
        fprintf('identifyer: %s\n',ERR.identifier);
        fprintf('message:    %s\n',ERR.message);
     end

    %% ===== THE WEL-file ====================================
    try
        i=strmatchi('WEL',nam(:,1));
        wel.NPER=dis.NPER;
        wel=readBCN([pth nam{i,3}],wel,'WEL');
    catch ERR
        fprintf('No WEL in %s\n',namfile);
        fprintf('identifyer: %s\n',ERR.identifier);
        fprintf('message:    %s\n',ERR.message);
     end

    %% ===== THE GHB-file ====================================
    try
        i=strmatchi('GHB',nam(:,1));
        ghb.NPER=dis.NPER;
        ghb=readBCN([pth nam{i,3}],ghb,'GHB');
    catch ERR
        fprintf('No GHB in %s\n',namfile);
        fprintf('identifyer: %s\n',ERR.identifier);
        fprintf('message:    %s\n',ERR.message);
    end

    %% ===== THE DRN-file ====================================
    try
        i=strmatchi('DRN',nam(:,1));
        drn.NPER=dis.NPER;
        drn=readBCN([pth nam{i,3}],drn,'DRN');
    catch ERR
        fprintf('No DRN in %s\n',namfile);
        fprintf('identifyer: %s\n',ERR.identifier);
        fprintf('message:    %s\n',ERR.message);
     end

    %% ===== THE RIV-file (river package) ===============================
    try
        i=strmatchi('RIV',nam(:,1));
        riv.NPER=dis.NPER;
        riv=readBCN([pth nam{i,3}],riv,'RIV');
    catch ERR
        fprintf('No RIV in %s\n',namfile);
        fprintf('identifyer: %s\n',ERR.identifier);
        fprintf('message:    %s\n',ERR.message);
    end

    %% ===== THE CHD-file FILE (constant head boundary package) =====
    try
        i=strmatchi('CHD',nam(:,1));
        chd.NPER=dis.NPER;
        chd=readBCN([pth nam{i,3}],chd,'CHD');
    catch ERR
        fprintf('No CHD in %s\n',namfile);
        fprintf('identifyer: %s\n',ERR.identifier);
        fprintf('message:    %s\n',ERR.message);
     end

    %% THE BTN-file (Basic Transprot Process for MT3D)
    try
        i=strmatchi('BTN',nam(:,1));
        btn=readBTN([pth nam{i,3}]); 
    catch ERR
        fprintf('No BTN in %s\n',namfile);
        fprintf('identifyer: %s\n',ERR.identifier);
        fprintf('message:    %s\n',ERR.message);
    end

    %% ===== THE ADV-file (Advection process) ==========
    try
        i=strmatchi('ADV',nam(:,1));
        adv=readADV([pth nam{i,3}]);
    catch ERR
        fprintf('No AVD in %s\n',namfile);
        fprintf('identifyer: %s\n',ERR.identifier);
        fprintf('message:    %s\n',ERR.message);
    end
                
    %% ===== THE DSP-file (Dispersion process package) =============
    try
        i=strmatchi('DSP',nam(:,1));
        dsp=readDSP([pth nam{i,3}]);
    catch ERR
        fprintf('No DSP in %s\n',namfile);
        fprintf('identifyer: %s\n',ERR.identifier);
        fprintf('message:    %s\n',ERR.message);
    end

    %% ===== THE SSM-file (Source-sink mixing process package) =====
    try
        ssm.NROW=btn.NROW; ssm.NCOL=btn.NCOL; ssm.NLAY=btn.NLAY; ssm.NPER=btn.NPER;
        i=strmatchi('SSM',nam(:,1));
        ssm=readSSM([pth nam{i,3}],btn);  % btn is necessary because of NROW NCOL NLAY and NPER
    catch ERR
        fprintf('No SSM in %s\n',namfile);
        fprintf('identifyer: %s\n',ERR.identifier);
        fprintf('message:    %s\n',ERR.message);
    end
    %% ===== THE rct-file (chemical reaction package) ==============
    try
        i=strmatchi('RCT',nam(:,1));
        rct=readRCT([pth nam{i,3}]);
    catch ERR
        fprintf('No RCT in %s\n',namfile);
        fprintf('identifyer: %s\n',ERR.identifier);
        fprintf('message:    %s\n',ERR.message);
    end

    %% ======================== SWI input ==========================
    try
        i=strmatchi('SWI',nam(:,1));
        swi.NROW=dis.NROW; swi.NCOL=dis.NCOL; swi.NLAY=dis.NLAY;
        swi=readSWI([pth nam{i,3}],swi);
    catch ERR
        fprintf('No SWI in %s\n',namfile);
        fprintf('identifyer: %s\n',ERR.identifier);
        fprintf('message:    %s\n',ERR.message);
    end
        

    %% ======================== CFP input ==========================
    try
        i=strmatchi('CFP',nam(:,1));
        cfp=readCFP([pth nam{i,3}]);
    catch ERR
        fprintf('No CFP in %s\n',namfile);
        fprintf('identifyer: %s\n',ERR.identifier);
        fprintf('message:    %s\n',ERR.message);
    end

    %% ======================== CRCH input ==========================
    try
        i=strmatchi('CRCH',nam(:,1));
        crch.NPER=dis.NPER; crch.NNODES=cfp.NNODES;
        crch=readCRCH([pth nam{i,3}],crch);
    catch ERR
        fprintf('No CRCH in %s\n',namfile);
        fprintf('identifyer: %s\n',ERR.identifier);
        fprintf('message:    %s\n',ERR.message);
    end

    %% ======================== CRCH input ==========================
    try
        i=strmatchi('COC',nam(:,1));
        coc=readCOC([pth nam{i,3}]);
    catch ERR
        fprintf('No COC in %s\n',namfile);
        fprintf('identifyer: %s\n',ERR.identifier);
        fprintf('message:    %s\n',ERR.message);
    end

end