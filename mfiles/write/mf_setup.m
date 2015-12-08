%function mf_setup
%MF_SETUP generates all input files for models (Backbone of mfLab)
%
% Example / USAGE:
%     mf_setup
%
%   This generic script generates inpt files for MODFLOW, MT3D and SEAWAT
%   and an increasing number of other packages for the MODFLOW suite.
%   First the paths are set to the directories containing the mfiles needed
%   by mflab.
%
%   In a second step mf_adapt in run locally, which is a script that
%   generates the actual model.
%   The matrices generated in mf_adapt are stored in basename.mat for
%   subsequent retrieval during the generation of the input files for the
%   groundwater models and as a backup.
%
%   The basename of the curent model is stored in basename.mat. Default is
%   the name of the local directory.
%   By default all generated filenames consist of basename plus an extension
%   denoting the function of the file.
%
%   mf_setup will use basename.xls, an excel file as a parameter container.
%   It is used to read the packages )to be invoked and the parameters that
%   need to be set in the simulations. These parameters are necessary for
%   the input files.
%
%   See the documentaion for more detailed information.
%
% TO 090101..130430
%
% Copyright 2009 2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

clear variables

%% Set paths, as mf_adapt may need them

fprintf('======= Start construction of your model ============\n');

fprintf('Your current HOME directory is\n%s\n',pwd);

% make a shortcut in the shortcut menu to set the path to the
% directories where the mfiles of mflab reside
% by clicking that shortcut after loading modflow
% Matlab is aware of these folders and the functions in them.
%     
% if ismac
%     MODELS='/Users/Theo/GRWMODELS/';
% else
%     MODELS='Z:\tolsthoorn On My Mac\GRWMODELS\';
% end
%
% MFLAB=[MODELS,'mflab' filesep]; % mflab home directory
% addpath([MFLAB 'write']);
% addpath([MFLAB 'read' ]);
% addpath([MFLAB 'gridcoords']);
% addpath([MFLAB 'etc']);
% addpath([MFLAB 'fdm']);

%% =====RUN MMF_ADAPT ===================================================
%  MF_ADAPT is a local m-file that may be used to the model matrices
%  created or to make a model completely
%  also in mf_adapt you set the paths for your local model

if exist('mf_adapt.m','file'),
    fprintf('Running model construction script mf_adapt..\n');
    mf_adapt;
elseif exist('mf_build.m','file');
    fprintf('Running model construction script mf_build..\n');
    mf_build;
else
    error(['mfLab:' mfilename ':no_mf_adapt_or_mf_build'],...
        '%s Can''t find mf_adapt.m nor mf_build.m\n',mfilename);
end

%% Save the model so we can generate input files

fprintf('Saving data to <<%s.mat>>\n',basename);

save('name.mat','basename');

S='save([basename,''.mat''],''IBOUND'',''STRTHD''';  % always

if exist('rGr','var'), xGr=rGr; end


%% Get the parameter types to grab wells, mnw, gridObj etc

clear ans;

%% Combine objects

mf_variables = whos;

Ipoint= strmatchi('pointObj',{mf_variables.class},'empty');
Iline = strmatchi('lineObj' ,{mf_variables.class},'empty');
Iarea = strmatchi('area2Obj',{mf_variables.class},'empty');

%% getting the user's gridObj if any and change its name to GRID,
%  save the users' grid along ith the other variables so it is not lost
%  when reloading basename both the user's grid and GRID will be present. 
myGridName='GRID';

clear GRID;

igrid=strmatchi('gridObj',{mf_variables.class},'exact'); % look for variables of class gridObj
if ~igrid  % if none, generate one using xGr,yGr,zGr and possibly LAYCBD
    if ~exist('xGr','var') || ~exist('yGr','var')
        error(['%s: If you did not define a grid by calling the gridObj, then\n',...
               'you must have xGr and yGr vectors in the workspace'],mfilename);
    end
    
    % There is no grid here. So this info must come from the workspace
    if ~exist('LAYCBD','var') || isempty(LAYCBD), LAYCBD=0; end
    if ~exist('AXIAL' ,'var') || isempty(AXIAL),  AXIAL=false; end
    if ~exist('MINDZ' ,'var') || isempty(MINDZ),  MINDZ=0.001; end
    
    try % zGr
        GRID=gridObj(xGr,yGr,zGr,LAYCBD,MINDZ,AXIAL);
    catch ME
        try % z
            GRID=gridObj(xGr,yGr,z,LAYCBD,MINDZ,AXIAL);
        catch ME
            try % Z
                GRID=gridObj(xGr,yGr,Z,LAYCBD,MINDZ,AXIAL);
            catch ME
                if ~exist('zGr','var') && ~exist('z','var') && ~exist('Z','var')
                    error(['%s: %s\nIf you did not define a grid by calling the gridObj, then\n',...
                        'you must define the model layer tops and bottoms throug zGr, z or Z in the workspace.'],...
                        mfilename,ME.message);
                end
            end
        end
    end
    
    yourGridName=myGridName;
    
elseif numel(igrid)==1 % if there is only one, use it
    
    yourGridName=mf_variables(igrid).name;
    
%     warning('on','mf_setup:GRID:useYourGridObj');
%     warning('mf_setup:GRID:useYourGridObj',...
%         ['mf_setup found %s class type variable: <<%s>>.\n',...
%         'mf_setup will use it for all grid requirements and grid info when generating the model input files\n.',...
%         'Make sure this is what you intended. Otherwise remove the gridObj at the end of mf_adapt.\n',...
%         'mfLab will temporarily rename it to <<%s>> and reconvert it to <<%s>> at the end of mf_setup.'],...
%         mf_variables(igrid).class,yourGridName,myGridName,yourGridName);
    
    eval(['GRID=' mf_variables(igrid).name]);

else
    s=mf_variables(igrid(1)).name;
    for i=2:numel(igrid), s=[s ', ' mf_variables(igrid(i)).name]; end  %#ok
    error('mf_setup:gridObj:tooManyGridObj',...
    ['you have <<%d>> variables of class ''gridObj'' defined in the workspace:\n',...
     '{ %s }\n',...
     'this is ambiguous, make sure only the correct one exists upon exiting mf_adapt.\n'],...
    numel(igrid),s);
end

save yourgridname yourGridName

% Hydraulic parameters for LPF package
if exist('LAYVKA','var'), S=[S,',''LAYVKA''']; LAYVKA = arrayCheck('LAYVKA',LAYVKA); end  % LPF kv/kh if LAYVKA ~=0
if exist('LAYCBD','var'), S=[S,',''LAYCBD''']; LAYCBD = arrayCheck('LAYCBD',LAYCBD); end  % for LPF
if exist('LAYCON','var'), S=[S,',''LAYCON''']; end  % for LPF
if exist('LAYTYP','var'), S=[S,',''LAYTYP''']; end  % for LPF
if exist('LAYAVG','var'), S=[S,',''LAYAVG''']; end  % for LPF
if exist('LAYWET','var'), S=[S,',''LAYWET''']; end  % for LPF
if exist('WETDRY','var'), S=[S,',''WETDRY''']; end  % for LPF

if exist('HK'    ,'var'), S=[S,',''HK'''    ]; HK    = arrayCheck('HK'   ,HK);    end  % LPF (alias for HY in BCF)
if exist('CHANI' ,'var'), S=[S,',''CHANI''' ]; CHANI = arrayCheck('CHANI',CHANI); end  % LPF ky/kx
if exist('HANI'  ,'var'), S=[S,',''HANI'''  ]; HANI  = arrayCheck('HANI' ,HANI);  end  % LPF ky/kx
if exist('VK'    ,'var'), S=[S,',''VK'''    ]; VK    = arrayCheck('VK'   ,VK  );  end  % LPF
if exist('VKA'   ,'var'), S=[S,',''VKA'''   ]; VKA   = arrayCheck('VKA'  ,VKA );  end  % LPF kv/kh if LAYVKA ~=0
if exist('SS'    ,'var'), S=[S,',''SS'''    ]; SS    = arrayCheck('SS'   ,SS  );  end  % for LPF
if exist('SY'    ,'var'), S=[S,',''SY'''    ]; SY    = arrayCheck('SY'   ,SY  );  end  % for LPF
if exist('VKCB'  ,'var'), S=[S,',''VKCB'''  ]; VKCB  = arrayCheck('VKCB' ,VKCB);  end  % for LPF

% Hydraulic parameters for BCF package
if exist('HY'    ,'var'), S=[S,',''HY'''    ]; HY    = arrayCheck('HY'   ,HY);    end  % for BCF=KH
if exist('TRAN'  ,'var'), S=[S,',''TRAN'''  ]; TRAN  = arrayCheck('TRAN' ,TRAN);  end  % for BCF
if exist('VCONT' ,'var'), S=[S,',''VCONT''' ]; VCONT = arrayCheck('VCONT',VCONT); end  % for BCF
if exist('SF1'   ,'var'), S=[S,',''SF1'''   ]; SF1   = arrayCheck('SF1'  ,SF1);   end  % for BCF
if exist('SF2'   ,'var'), S=[S,',''SF2'''   ]; SF2   = arrayCheck('SF2'  ,SF2);   end  % for BCF

if exist('HFB'   ,'var'); S=[S,',''HFB'''  ];  end  % horizontal flow barrier

% Recharge and evaptranspiration packages
if exist('RECH'   ,'var'), S=[S,',''RECH''' ]; end
if exist('IRCH'   ,'var'), S=[S,',''IRCH''' ]; end
if exist('SURF'   ,'var'), S=[S,',''SURF''' ]; end
if exist('EVTR'   ,'var'), S=[S,',''EVTR''' ]; end
if exist('EXDP'   ,'var'), S=[S,',''EXDP''' ]; end
if exist('IEVT'   ,'var'), S=[S,',''IEVT''' ]; end

% ETS evapotranspiration package ( multi segment )
if exist('PXDP'   ,'var'), S=[S,',''PXDP''' ]; end
if exist('PETM'   ,'var'), S=[S,',''PETM''' ]; end

% Boundary conditions for MODFLOW
if exist('WEL'   ,'var'), S=[S,',''WEL'''   ]; end
if exist('MNW'   ,'var'), S=[S,',''MNW'''   ]; end
if exist('DRN'   ,'var'), S=[S,',''DRN'''   ]; end
if exist('RIV'   ,'var'), S=[S,',''RIV'''   ]; end
if exist('GHB'   ,'var'), S=[S,',''GHB'''   ]; end
if exist('CHD'   ,'var'), S=[S,',''CHD'''   ]; end

% save pointObj, lineObj, area2Obj
for i=[Ipoint Iline Iarea]
    S=[S,',''' mf_variables(i).name, '''']; %#ok
end
    
%% Get all data of type wellObj

iWell = strmatchi('wellObj',{mf_variables.class},'empty');
for iw = iWell(:)', S=[S,',''',mf_variables(iw).name,'''']; end %#ok

iWell = strmatchi('MNW1Obj',{mf_variables.class},'empty');
for iw=iWell
    if strcmp('mnw1',mf_variables(iw).name)
        error('%s: the name <<%s>> is not allowed, it interferes with mf_setup',...
            mfilename,mf_variables(iw).name);
    end
    S=[S,',''',mf_variables(iw).name,'''']; %#ok
end

iWell = strmatchi('MNW2Obj',{mf_variables.class},'empty');
for iw=iWell
    if strcmp('mnw2',mf_variables(iw).name)
        error('%s: the name <<%s>> is not allowed, it interferes with mf_setup',...
            mfilename,mf_variables(iw).name);
    end
    S=[S,',''',mf_variables(iw).name,'''']; %#ok
end

% MODPATH, add particle groups irrespective of their names
I = strmatchi('mpath_particleGroupObj',{mf_variables.class},'empty');
for i=I
    fprintf('mpath_particleGroupObj = %s\n',mf_variables(i).name);
end
if numel(I)>1
    error(['%s: the number of mpath_particleGroupObj is <<%d>>.',...
    'This is ambiguous, join them to a single array of such objects in',...
    'mf_adapt (or mf_build)'],mfilename,numel(I));
end
if ~isempty(I)
    particleGroupNm = mf_variables(I).name;
    S=[S,',''particleGroupNm'',''' particleGroupNm ''''];
end

if exist('ZONE',         'var'), S=[S,',''ZONE'''         ]; ZONE          = arrayCheck('ZONE',ZONE); end  % MODPATH stop zone array
% 'stopZone' --> in worksheet MPATH, not in workSpace
if exist('RETARDATION'  ,'var'), S=[S,',''RETARDATION'''  ]; RETARDATION   = arrayCheck('RETARDATION'  ,RETARDATION);   end  % MODPATH retardation array
if exist('RETARDATIONCB','var'), S=[S,',''RETARDATIONCB''']; RETARDATIONCB = arrayCheck('RETARDATIONCB',RETARDATIONCB); end  % MODPATH same for conf. beds
if exist('BUDGETCELLS','var'),   S=[S,',''BUDGETCELLS'''  ]; end  % MODPATH cell Idx for wich budget is checked

% Boundary conditions transport MT3DMS and Seawat
if exist('ICBUND','var'), S=[S,',''ICBUND''']; ICBUND = arrayCheck('ICBUND',ICBUND); end  % MT3D
if exist('STCONC','var'), S=[S,',''STCONC''']; end  % MT3D BTN
if exist('PNTSRC','var'), S=[S,',''PNTSRC''']; end  % MT3D SSM
if exist('BTNOBS','var'), S=[S,',''BTNOBS''']; end  % MT3D BTN

% Recharge and EVT concentrations in SSM
if exist('CRCH','var'), S=[S,',''CRCH''']; end  % MT3D SSM
if exist('CEVT','var'), S=[S,',''CEVT''']; end  % MT3D SSM

% Necessary for transport, effective porosity
% if two alternatives exist for this variable mfLab preferes PEFF
if exist('POROSITY','var') && ~exist('PEFF','var'), PEFF = arrayCheck('POROSITY',POROSITY); clear POROSITY; end
if exist('PRSITY','var')   && ~exist('PEFF','var'), PEFF = arrayCheck('PRSITY'  ,PRSITY);   clear PRSITY;   end
if exist('POR'   ,'var')   && ~exist('PEFF','var'), PEFF = arrayCheck('POR'     ,POR);      clear POR;      end

if exist('PEFF'  ,'var'),  S=[S,',''PEFF'''  ]; PEFF   = arrayCheck('PEFF'  ,PEFF);   end  % MT3D BTN
if exist('PORCB'  ,'var'), S=[S,',''PORCB'''];  PORCB  = arrayCheck('PORCB' ,PORCB);  end  % PEFF of confining bed
if exist('DMCOEF','var'),  S=[S,',''DMCOEF''']; DMCOEF = arrayCheck('DMCOEF',DMCOEF); end  % MT3D DSP
if exist('AL',    'var'),  S=[S,',''AL'''    ];  end  % MT3D DSP

if exist('RHOB'   ,'var'), S=[S,',''RHOB'''   ]; RHOB   = arrayCheck('RHOB'   ,RHOB);    end % rct
if exist('PRSITY2','var'), S=[S,',''PRSITY2''']; PRSITY2= arrayCheck('PRSITY2',PRSITY2); end % rct
if exist('SRCONC' ,'var'), S=[S,',''SRCONC''' ]; SRCONC = arrayCheck('SRCONC' ,SRCONC);  end % rct
if exist('SP1'    ,'var'), S=[S,',''SP1'''    ]; SP1    = arrayCheck('SP1'    ,SP1);     end % rct 
if exist('SP2'    ,'var'), S=[S,',''SP2'''    ]; SP2    = arrayCheck('SP2'    ,SP2);     end % rct
if exist('RC1'    ,'var'), S=[S,',''RC1'''    ]; RC1    = arrayCheck('RC1'    ,RC1);     end % rct
if exist('RC2'    ,'var'), S=[S,',''RC2'''    ]; RC2    = arrayCheck('RC2'    ,RC2);     end % rct

% SEAWAT (direct reading in of DENSE and VISC probably never going to be used
if exist('DENSE' ,'var'), S=[S,',''DENSE''' ]; DENSE = arrayCheck('DENSE',DENSE);  end  % SEAWAT
if exist('VISC'  ,'var'), S=[S,',''VISC'''  ]; VISC  = arrayCheck('VISC', VISC);   end  % SEAWAT

% Additional options for CHD in relation to Seawat, see Langevin e.a. 2008,
% p12-15 & 22
if exist('CHDDENSOPT','var'), S=[S,',''CHDDENSOPT''']; end % CHD-SEAWAT
if exist('CHDDENS','var'), error('Use CHDDEN instead of CHDDENS as described in the SEAWAT manual.'); end
if exist('CHDDEN','var')   , S=[S,',''CHDDEN'''   ]; end % CHD-SEAWAT

% Necessary for SWI, the Saltwater Instrusion Package
% p12-15 & 22
if exist('ISOURCE','var'), S=[S,',''ISOURCE''']; ISOURCE = arrayCheck('ISOURCE',ISOURCE); end % SWI
if exist('ZETA'   ,'var'), S=[S,',''ZETA'''   ]; end % SWI
if exist('SSZ'    ,'var'), S=[S,',''SSZ'''    ]; SSZ     = arrayCheck('SSZ'    ,SSZ);     end % SWI

% To allow execuring specific mfiles after mf_setup
if exist('AFTERMFSETUP','var'), S=[S,',''AFTERMFSETUP''']; end % PEST

% Axisymmetric flow
if     exist('RADIAL','var'),       AXIAL=RADIAL;        clear RADIAL;
elseif exist('AXISYMMETRIC','var'), AXIAL=AXISYMMETRIC;  clear AXISYMMETRIC;   
elseif exist('AXIAL','var');        % skip
else   AXIAL = false;
end

AXIAL = AXIAL || GRID.AXIAL;

% Set AXIAL in grid and make sure all internal grid procedures are honored
if GRID.AXIAL
    % all DY are 1 if AXIAL !! This is to protect the user from errors in
    % case DY~=1, in which case MODLFOW will compute wrong recharge etc.
    % It must use DY==1 internally to be consistent.
    GRID=gridObj(GRID.xGr, GRID.yGr(1)+(0:GRID.Ny), GRID.Z, GRID.LAYCBD, GRID.MINDZ, AXIAL);
else
    GRID=gridObj(GRID.xGr, GRID.yGr, GRID.Z, GRID.LAYCBD, GRID.MINDZ, AXIAL);
end
if exist('RUNBACKGROUND','var'), BACKGROUND=RUNBACKGROUND; end

if ~exist('BACKGROUND','var'), BACKGROUND=0; end
if ~exist('GREP','var')      , GREP=''     ; end

S=[S,',''AXIAL'',''BACKGROUND'',''GREP'',''',myGridName,''',''',yourGridName,''''];

S=[S,');'];

try
    eval(S); fprintf('%s\n',S);
catch ME
    fprintf(['\n\nTrouble saving the variables from the Matlab workspace.\n',...
        'One of the required variables seems to be missing in the Matlab workspace.\n',...
        'This variable is necessary for the correct working of MODFLOW\n',...
        'or one of its related programs. Check if the variable has been provided.\n',...
        'Also check its correct spelling. See the error message that follows\n',...
        'for the name of the variable that is missing:\n']);
    throw(ME)
end

fprintf('.... Saved !\n');

%%  =====START GENERATING MODFLOW etc INPUT FILES ============================
% Cleaning up the workspace removes cluttering variables and is important
% to make sure later on that you are working with the correct variables
% saved in <<basename.mat>> and that you saved all variables necessary
% to run the model

fprintf('Cleanup workspace ....');

clear variables                % clean workspace, but leave debug settings
%close all                      % close all figures
load('name'); load([basename '.mat']);  % load data computed in mf_adapt and just saved

fprintf('... workspace cleaned and <<%s.mat>> reloaded.\n',basename);

%% ======= IF NECESSRY ADAPT TO YOUR OWN SITUATION ====================

tic;   % tick is to measure elapsed time from this point toc to present it

% Path to the executables. I put them all in MODELS/bin
% because I have several differntly compiled versions. It's your
% choise, anyway set the parameters MODFLOW MT3D etc down below to
% their actual lcoations on your hard drive.

% NOTICE the " " in the paths below to manage spaces in file names
% in system command used to launch the external executable later on

%% Get workspace variables

clear ans
mf_variables = whos; % variables from Matlab's workspace

fprintf('Defining paths to your excecutables\n');

% Legal models, the list may be extended

LEGAL_MODELS={'MF2000','MF2005','MF2007','MF2005NWT','MF2KASP','MFSWI','SEAWAT','MT3MD','MODPATH'};

setExecutables;

%% STARTING WITH THE ACTUAL DATA FOR THIS MODEL =============

% ===== The EXCEL workbook with all parameter settings ===============
fprintf('Basename current model is ''%s''\n',basename);

if     exist([basename,'.xls' ],'file'), XLSF=[basename,'.xls'];
elseif exist([basename,'.xlsx'],'file'), XLSF=[basename,'.xlsx'];
else
    error('Can''t find accompanying Excel workbook <<%s>> nor <<%s>>',...
    [basename,'.xls'],[basename,'.xlsx']);
end
fprintf('Getting nam file data from %s\n',XLSF);

%% xlsread mac is different from xlsread windows !!!

if ismac
    warning('off','MATLAB:xlsread:Mode')
    try
        [Nnum,Ntxt]=xlsread(XLSF,'NAM','','basic');
    catch ME
        if ~isempty(strfind(ME.identifier,'WorksheetNotFound'))
                error(['%s ''NAM''\n',...
                    'Save the worksheet as an Excel 5.0/95 Workbook (.xls)',...
                    ' so that Matlab can read it on non-Windows systems\n'],...
                    ME.message);
        end
    end
    warning('on','MATLAB:xlsread:Mode')
else
    [Nnum,Ntxt]=xlsread(XLSF,'NAM');    
end

% Get the data for the *.nam file from the NAM worksheet in the XSLF workbook 

%% Another incompatibility between xlsread mac and windows
if ismac   %  size(Nnum,2)==size(Ntxt,2)
    Nnum=Nnum(2:end,:);  % cut off header line
    Ntxt=Ntxt(2:end,:);  % cut off header line
else
    % on pc Nnum needs not be cut contrary to Mac
    if size(Nnum,1)>1 && all(isnan(Nnum(1,:))), Nnum=Nnum(2:end,:); end
    Ntxt=Ntxt(2:end,:);  % cut off header line    
end

%%
nam.PCKG=Ntxt(:,1);     % List of packages
nam.EXT =Ntxt(:,3);     % corresponding file extenstions
nam.UNIT=Nnum(:,1);     % corresponding file units
nam.SCEN=Nnum(:,3);
		% =0 package not inlcuded
		% >0 value is scenario number. This corresponds with data column
	    %     in the worksheets MFLOW, MT3D and SWT next to the package/parnam column
	    %     this column must be present when package is processed!!
        % <0  for VDF run seawat but do not run package
        % <0  for SWI run swi    but do not run package

% ====UNIT numbers for the output files (hds, ddn, bgt,zta) ==============
% These units are associate with the given extension "hds" "ddn" "bgt" "zta"
% and not with the specified package names

%% Sort out which model to use, default MF2K
%  This is done on the NAM worksheet using pseudo packages MF2, MF5 and
%  ASP, SWT, MPTH

nam.MODEL   = {'MF2000'};
nam.mdlpath = {MF2000};

mflowWasRunByUser = false;

i=strmatchi('MF2000',nam.PCKG,'exact');
if (i>0 && nam.SCEN(i)>0),
    fprintf('NAM worksheet specifies use of MF2000\n');
    nam.MODEL{end}   = 'MF2000';
    nam.mdlpath{end,1} =  MF2000;
    mflowWasRunByUser = true;
end

% overrule if set:
i=strmatchi('MF2005',nam.PCKG,'exact');
if (i>0 && nam.SCEN(i)>0),
    fprintf('NAM worksheet specifies use of MF2005 instead of MF2k\n');
    nam.MODEL{end}   = 'MF2005';
    nam.mdlpath{end,1} =  MF2005;
    mflowWasRunByUser = true;
end

% overrule if set
i=strmatchi('MF2007',nam.PCKG,'exact');
if (i>0 && nam.SCEN(i)>0),
    fprintf('NAM worksheet specifies use of MF2007 (Conduit Flow Package, CFP) instead of MF2k\n');
    nam.MODEL{end}   = 'MF2007';
    nam.mdlpath{end,1} =  MF2007;
    nam.PCKG(i)=[]; nam.EXT(i)=[]; nam.UNIT(i)=[]; nam.SCEN(i)=[];
    mflowWasRunByUser = true;
end

% overrule if set
i=strmatchi('MF2005NWT',nam.PCKG,'exact');
if (i>0 && nam.SCEN(i)>0),
    fprintf('NAM worksheet specifies use of MF2005NTW (2005 with Newton Iteragin Proces) instead of MF2k\n');
    nam.MODEL{end}   = 'MF2005NWT';
    nam.mdlpath{end,1} =  MF2005NWT;
    nam.PCKG(i)=[]; nam.EXT(i)=[]; nam.UNIT(i)=[]; nam.SCEN(i)=[];
    mflowWasRunByUser = true;
end

% overrule if set
i1 = strmatchi('MF2KASP',nam.PCKG,'exact');
if (i1>0 && nam.SCEN(i1)>0)
    if ispc
        nam.MODEL{end}   = 'MF2KASP';
        nam.mdlpath{end,1} =  MF2ASP;
        warning('on','mf_setup:modelChoice:ambiguous');
        warning('mf_setup:modelChoice:ambiguous',...
            'Will run MF2ASP, if not desired, switch off model MF2ASP and package ASP!\n');        
    else
        error('mf_setup:modelChoice:incopatible',...
            ['MF2KASP specified in nam sheet, but John Doherty''s MF2KASP is only available under windows.\n',...
            'choose another model to run (see list in top of NAM sheet)!']);
    end
    fprintf('NAM worksheet specifies use of MF2ASP model instead of MF2K\n');
    mflowWasRunByUser = true;
end

% possibly overrule if set
% SWI is chosen if:
%     SWI pacakge is on and no model has as yet been selected
%  or SWI package is on and MF2000 has been selected.
i1=strmatchi('MFSWI',nam.PCKG,'exact');
i2=strmatchi('SWI'  ,nam.PCKG,'exact');
if (i1>0 && nam.SCEN(i1)>0) || ...
   (i2>0 && nam.SCEN(i2)>0 && strmatchi('MF2000',nam.MODEL)) &&~strmatchi('MT3DMS',nam.MODEL)
        nam.MODEL{end}     ='MFSWI';
        nam.mdlpath{end,1}   = MFSWI;
    mflowWasRunByUser = true;
end

% possibly overrule by SEAWAT if
%    SEAWAT is on
% or ~MODEL & VDF
i1=strmatchi('SEAWAT',nam.PCKG,'exact');
i2=strmatchi('BTN',nam.PCKG,'exact');
i3=strmatchi('VDF',nam.PCKG,'exact');
if  (i1>0 && nam.SCEN(i1)>0) || ...
    isempty(nam.MODEL) && (i2>0 && nam.SCEN(i2)>0) && (i3>0 && abs(nam.SCEN(i3))>0)
    fprintf('NAM worksheet specifies use of SEAWAT\n');
    nam.MODEL{end}   ='SEAWAT';
    nam.mdlpath{end,1} = SEAWAT;
end

% possibly overrule by MT3DMS if
%      MT3DMS is on
% or   ~MODEL && (LMT6 && BTN are on)
% and
% We have to run MODFLOW before running MT3DMS
i1=strmatchi('MT3DMS',nam.PCKG);
i2=strmatchi('BTN'   ,nam.PCKG);
i3=strmatchi('LMT6'  ,nam.PCKG);
if  (i1>0 && nam.SCEN(i1)>0) || ...
    (isempty(nam.MODEL) && (i2>0 && nam.SCEN(i2)>0) && (i3>0 && nam.SCEN(i3)>0))
    fprintf('NAM worksheet specifies use of MT3DMS instead of MF2k\n');
    if ~mflowWasRunByUser
        nam.MODEL{end}       = 'MT3DMS';
        nam.mdlpath{end,1}   =  MT3DMS;
    else
        nam.MODEL{end+1}       = 'MT3DMS';
        nam.mdlpath{end+1,1}   =  MT3DMS;
    end
end
i=strmatchi('MODPATH',nam.PCKG,'exact');
if i>0 && nam.SCEN(i)>0
    fprintf('NAM worksheet specifies generating intput for MODPATH\n');
    nam.MODEL{end+1}   ='MODPATH';
    nam.mdlpath{end+1,1} = MODPATH;
end

%% Check to see that we can only have one MF model

fprintf('You selected the following MODEL residing on the following path\n');
for iModel=1:numel(nam.MODEL)
    fprintf('Model: %s\n',nam.MODEL{iModel});
    fprintf('Path:  %s\n',nam.mdlpath{iModel});
end
        
%% Start generation of input files

fprintf('\nStarting generation of input files ....\n');

%% Legal packages, the list may be extented

% to check below whether a legal solver is on:
solvers = {'PCG','DE4','SIP','SOR','NWT','PGCN'};

% Specific for MODFLOW
modflowPackages={...
    'BAS6','DIS',...
    'BCF6','LPF',...
    'RCH','EVT','ETS',...
    'WEL','MNW1','MNW2','RIV','GHB','DRN','CHD',...
    'HFB6',...
    'PES','SEN','OBS','HOB','LMT6','UMT',...
    'CFP','COC','CRCH',...
    'PCG','DE4','SIP','SOR','PCGN'...
    'OC','DATA ','DATA(BINARY)',...
    };

transportpackages={'BTN','ADV','DSP','SSM','GCG','RCT'};

% modpath version 6 only (2012)
mpath_packages={'MPBAS','DIS','HEAD','BUDGET','DATA','DATA(BINARY)'};
        
nam.LegalPCKG=cell(numel(nam.MODEL),1);

for iModel=1:numel(nam.MODEL)
    switch nam.MODEL{iModel}
        case 'MF2000',
            % 'GLOBAL' removed TO 121210, obsolete
            nam.LegalPCKG{iModel} = horzcat('LIST',modflowPackages);
        case 'MF2005',
            nam.LegalPCKG{iModel} = horzcat('LIST',modflowPackages,'SWI2');
        case 'MF2007',
            nam.LegalPCKG{iModel} = horzcat('LIST',modflowPackages);
        case 'MF2005NWT'
            nam.LegalPCKG{iModel} = horzcat('LIST','NWT','UPW',modflowPackages);
        case 'MF2KASP'
            nam.LegalPCKG{iModel} = horzcat('LIST','ASP',modflowPackages);
        case 'MT3DMS'
            % if MODFLOW (any version) is to be run before MT3DMS, then make sure the
            % packages LMT6 and FTL are on for MF2000
            % this is regardless of the user's setting in the NAM
            % worksheet.
            imf = strmatchi('MF',nam.MODEL); % MF2000, MF2005 etc all starting with MF
            if imf,
                nam.LegalPCKG{imf}=horzcat('LMT6',nam.LegalPCKG{imf});
            end
            nam.LegalPCKG{iModel} = horzcat('LIST','FTL',transportpackages);
        case 'SEAWAT'
             nam.LegalPCKG{iModel} = horzcat('LIST','VDF','VSB',modflowPackages,transportpackages,'VDF','VSC');
        case 'MFSWI'
             nam.LegalPCKG{iModel} = horzcat('LIST',modflowPackages,'SWI');
        case 'MODPATH'
             nam.LegalPCKG{iModel} = horzcat(mpath_packages);
    end
end

% Make sure packages are unique, this must be done separately because of
% the measures taken for the MODFLOW-MT3DMS combination.
for iModel = numel(nam.MODEL)
    nam.LegalPCKG{iModel} = unique(nam.LegalPCKG{iModel});
end

%% REMOVE all lines with UNIT non numerical

ison= ~isnan(nam.UNIT) & nam.SCEN~=0;

if all(~ison), error('mf_setup: no packages selected'); end

nam.PCKG = nam.PCKG(ison);
nam.EXT  = nam.EXT (ison);
nam.UNIT = nam.UNIT(ison);
nam.SCEN = nam.SCEN(ison);

%% Check if unit numbers are unique

for i=1:numel(nam.UNIT)-1;
    if any(nam.UNIT(i)==nam.UNIT(i+1:end))
        error('mf_Setup:UNITnr:NotUnique',...
            'Unit number <<%d>> in NAM file is not unique !',nam.UNIT(i));
    end
end

%% Check uniqueness of file extensions

for i=1:numel(nam.EXT)-1
    if any(strmatchi(nam.EXT{i},nam.EXT(i+1:end)))~=0,
            error('File extension <<%s>> in nam file is not unique !',nam.EXT{i});
    end
end

%% switch off all illegal packages

allLegal = unique(horzcat(nam.LegalPCKG{:}));
islegal = zeros(size(nam.PCKG));
for i=1:numel(nam.PCKG),
    islegal(i)=strmatchi(nam.PCKG{i},allLegal,'exact');
end

nam.PCKG = nam.PCKG(islegal>0);
nam.EXT  = nam.EXT (islegal>0);
nam.UNIT = nam.UNIT(islegal>0);
nam.SCEN = nam.SCEN(islegal>0);

%% NAM Namefile INFO
fprintf('Preparing name file struct.\n');

% GENERATE THE NAM+BAT FILES FOR MODFLOW, MT3D AND SEAWAT
% output nam has some implicit packages added in case of MODPATH is used
% namely MAIN (and LOC)
nam = writeNAM(basename,nam);

%% See if any isnan(IBOUND), if so set corresponding 3D values equal to 0
% This is to prevent NaN in data passed to model input, which it cannot handle
IBOUND(isnan(IBOUND))=0;

%% SIMULATION INFO MODFLOW (always needed)
fprintf('Getting simulation parameters from worksheet %-8s in %s\n','MFLOW',XLSF);
[MFLOWparnams,MFLOWparvals,MFLOWtxtnams,MFLOWtxtvals]=getExcelData(XLSF,'MFLOW','Vertical');

HDRY   = MFLOWparvals(strmatchi('HDRY'  ,MFLOWparnams));
HNOFLO = MFLOWparvals(strmatchi('HNOFLO',MFLOWparnams));
save HINACT HDRY HNOFLO

%% SIMULATION INFO MT3D (if  MT3D is to be run)
if strmatchi('BTN',nam.PCKG) || strmatchi('MT3',nam.PCKG)
    fprintf('Getting simulation parameters from worksheet %-8s in %s\n','MT3D',XLSF);
    [MT3Dnams,MT3Dvals]=getExcelData(XLSF,'MT3D','Vertical');
    CINACT = MT3Dvals(strmatchi('CINACT',MT3Dnams));
    save CINACT CINACT;
end

%% SIMULATION INFO SEAWAT (if  SEAWAT is to be run)
if any(strmatchi({'VDF','VSC','SWT'},nam.PCKG))
    fprintf('Getting simulation parameters from worksheet %-8s in %s\n','SEAWAT',XLSF);
    [SEAWATparnams,SEAWATparvals]=getExcelData(XLSF,'SEAWAT','Vertical');
end

%% STRESS PERIOD INFO (always needed)
fprintf('Getting stress period data from    worksheet %-8s in %s\n','PER',XLSF);
[PERnams,PERvals,NPER]=getPeriods(XLSF);

if size(PERnams,2) ~= size(PERvals,2)
    error('mf_setup:getPeriods:NrOfHdrsDoesNotMatchColumns',...
        ['Number of headers does not match number of columns in PER sheet.\n',...
         'Verify that there are no empty columns or empty header labels,\n',...
         'and that there are no fommulas in the header labels.']);
end

%% LAYER INFO (always needed)

[LAYparnams,LAYparvals]=getLayers(XLSF,GRID);

if size(LAYparnams,2) ~= size(LAYparvals,2)
    error('mf_setup:getLayers:NrOfHdrsDoesNotMatchColumns',...
        ['Number of headers does not match number of columns in LAY sheet.\n',...
         'Verify that there are no empty columns or empty head labels,\n',...
         'and that there are no fommulas in the header labels.']);
end

if exist('LAYCBD','var')
    warning('on','mf_setup:LAYCBD:LgnoreSpreadSheet');
    warning('mf_setup:LAYCBD:LgnoreSpreadSheet',...
        'mf_setup: LAYCBD found in workspace; LAYCBD in worksheet LAY is ignored.\n');
end

%% PEST INFO (only needed with ASP and PES package)
if  any(strcmp('ASP',nam.PCKG)) && any(strcmp('PES',nam.PCKG))
    fprintf('Getting PEST   info from            worksheet %-8s in %s\n','PEST',XLSF);
    [PESTparnams,PESTparvals]=getExcelData(XLSF,'PEST','Vertical');
end

%% Check if the necessary DATA(BINARY) and DATA files are defined for hds
% and ddn (heads and drawdown)

ioc=strmatchi('OC',nam.PCKG,'exact');

I___UN={'HED' ,'DDN'};   % see I*UN in MFLOW worksheet
SV    ={'HDSV','DDSV'};  % see PER worksheet
for i=1:numel(I___UN)
    % do we want output for HEAD or DDN?
    
    if ioc>0 && nam.UNIT(ioc)>0 && any(PERvals(:,strmatchi(SV{i},PERnams)))
        
        % Then look for the corresponding output file in the NAM fle/worksheet
        % Is any DATA or DATA(BINARY) output file defined in NAM ???
        Idata=strmatchi('DATA',nam.PCKG);
        if ~Idata   % if none, throw error message and stop:
            error(['To write output as desired by the variable <<%s>> in the PER worksheet\n',...
                'you must specfy an output file for it in the NAM worksheet\n',...
                '(i.e. one line starting with DATA or DATA(BINARY))\n',...
                'which must be in accordance with the unit number for I%sUN in the\n',...
                'MFLOW worksheet. You may omit all output by switching off\n',...
                'the output package (OC) in the NAM file/worksheet.'],...
                SV{i},I___UN{i});
        end

        % What is the desired output unit specified in MOFLOW worksheet?
        iunit=MFLOWparvals(strmatchi(['I',I___UN{i},'UN'],MFLOWparnams,1));
        
        % is this the unit in one of the active data files in the NAM file/worksheet?
        if ~any(iunit==nam.UNIT(Idata))
            error(['%s: unit <<%d>> to output for <<I%sUN>> is not present in nam sheet or not on.\n',...
                 'REMEDY: Check the NAM sheet and MFLOW sheet for this unit number.\n',...
                 'Perhaps the binary data file for drawdown or similar is not switched on.\n'],...
                mfilename,iunit,I___UN{i});
        end
    end
end

%% Check if the necessary DATA(BINARY) and DATA files are defined

% These are the modules that may output cell-by-cell flow terms
% The concerned parameter names are obtained by preceding them with
% capital 'I' and appending 'CB'.
ICB={'LPF','BCF','RIV','DRN','GHB','CHD','WEL','MNW','EVT','ETS','RCH'};

% We will throw this error string in case an error is found

ioc=strmatchi('OC',nam.PCKG,'exact');

% if any output of cell-by-cell flow terms is desired (see ICBCFL in PER worksheet) 
if (ioc>0 && nam.SCEN(ioc)) && any(PERvals(:,strmatchi('ICBCFL',PERnams)))
    % Then we must find the desired unit number in one of the NAM-file lines
    % specifying output with DATA or DATA(BINARY) in the first column.

    Idata=strmatchi('DATA',nam.PCKG);  % find these DATA, DATA(BINARY) lines
    if ~Idata
        error(['To write cell-by-cell flow terms, you must specfy output files',...
            '(these are lines starting with DATA or DATA(BINARY))\n',...
            'in the NAM file/NAM worksheet in accordance with the unit numbers',...
            'for the packages I*CB in the MFLOW worksheet\n',...
            'If you do not want any cell by cell output, set ICBCFL=0 in PER worksheet.']);
    end
    
    % check for the individual packages that may write cell by cell output   
    for i=1:numel(ICB)  % all possible such packages
        % get the desired output unit number from MFLOW worksheet
        idx=strmatchi(['I',ICB{i},'CB'],MFLOWparnams,1);
        if idx==0, icb_unit=0; else icb_unit=MFLOWparvals(idx); end
        if icb_unit>0
            if ~any(icb_unit==nam.UNIT(Idata))
                error(['Check the NAM data: Missing DATA output unit number <<%d>>\n',...
                    'to write cell-by-cell flow terms for <<%s>> as specified in\n',...
                    'worksheet MFLOW !'],icb_unit,ICB{i});
            end
        end
    end
end

%% Start producing the model data

MXSS = 0; % counter max number of sources and sinks necesary for MT3DMS and SEAWAT

FREE = MFLOWparvals(strmatchi('FREE',MFLOWparnams),1)~=0;
if isempty(FREE), FREE=false; else FREE = logical(FREE); end

%% ===== THE BAS FILE ====================
bas.SCEN=namSCEN('BAS6',nam.PCKG,nam.SCEN);

if bas.SCEN
    fprintf('Generating basic struct\n');
    bas.unit=nam.UNIT(strmatchi('BAS6',nam.PCKG));
    bas.ext =nam.EXT {strmatchi('BAS6',nam.PCKG)};   
%    bas.NLAY=size(IBOUND,3); 
    bas.NLAY=GRID.Nlay;
    bas.HNOFLO =MFLOWparvals(strmatchi('HNOFLO',MFLOWparnams),bas.SCEN);
    bas.IBOUND = IBOUND;
    bas.STRTHD = STRTHD;  % initial heads zero
    bas.FREE   = FREE;
    bas.CHTOCH = MFLOWparvals(strmatchi('CHTOCH',MFLOWparnams),bas.SCEN)~=0;
    writeBAS6(basename,bas);
end

%% ===== THE DIS-file struct ==============
% STRESS PERIODS
dis.SCEN=namSCEN('DIS',nam.PCKG,nam.SCEN);

if dis.SCEN
    fprintf('Generating discretization struct\n');
    dis.unit=nam.UNIT(strmatchi('DIS',nam.PCKG));
    dis.ext =nam.EXT {strmatchi('DIS',nam.PCKG)};
    
    dis.GRID=GRID;  % not yet used in writedis (April 2012)
    dis.FREE=bas.FREE;
    
    dis.NPER=NPER;
        
    dis.ITMUNI =MFLOWparvals(strmatchi('ITMUNI',MFLOWparnams),dis.SCEN);
    dis.LENUNI =MFLOWparvals(strmatchi('LENUNI',MFLOWparnams),dis.SCEN);  % confining bed below this layer

    % Generating stress period info for discretizaton file: 
    dis.PERLEN = PERvals(:,strmatchi('PERLEN',   PERnams));
    dis.NSTP   = PERvals(:,strmatchi('NSTP',     PERnams));
    dis.TSMULT = PERvals(:,strmatchi('TSMULT',   PERnams));
    dis.isTran = PERvals(:,strmatchi('Transient',PERnams));
    
    writeDIS(basename,dis);
end

%% ===== THE HFB (Horizontal Barrier Package) file =============
hfb.SCEN=namSCEN('HFB',nam.PCKG,nam.SCEN);

if hfb.SCEN
    fprintf('Generating Horizontal Flow Barrier (hfb)struct\n');
   
    if ~exist('HFB','var')
        error(['variable HFB not defined in the workspace\n',...
               'It must be [iPer L R C R C Hydcdr]']);
    end
   
    hfb.FREE=bas.FREE;
    hfb.NPER=NPER;
    hfb.unit=nam.UNIT(strmatchi('HFB',nam.PCKG));
    hfb.ext =nam.EXT {strmatchi('HFB',nam.PCKG)};

    hfb.NLAY=GRID.Nlay;
    hfb.NROW=GRID.Ny;
    hfb.NCOL=GRID.Nx;

    %1
    hfb.NPHFB=0; % number of parameters, don't use parameters
    hfb.MXFB=0;  % max nr of HFB barriers to be defined through parameters
    hfb.NHFBNP=size(HFB,1);  % number of HBB barriers not defined through parameters
    %2
    %3
    %4
    hfb.HFB=HFB;
    if GRID.AXIAL
        if size(HFB,2)==6,
            R=0.5*abs(GRID.xm(HFB(:,3))+GRID.xm(HFB(:,5)));
        else
            R=0.5*abs(GRID.xm(HFB(:,4))+GRID.xm(HFB(:,7)));
        end
        hfb.HFB(:,end)=2*pi*R.*HFB(:,end);
    end
    %5
    %5
    hfb.NACTHFB=0; % number of active HFB parameters 
    writeHFB(basename,hfb);
end

%% ===== THE BCF-file =====================
bcf.SCEN=namSCEN('BCF6',nam.PCKG,nam.SCEN);

if bcf.SCEN    
    fprintf('Generating BCF struct and file\n');
    bcf.unit  =nam.UNIT(strmatchi('BCF6',nam.PCKG));
    bcf.ext   =nam.EXT {strmatchi('BCF6',nam.PCKG)};
    
    bcf.FREE = bas.FREE;
    
    bcf.GRID = GRID;
    
    bcf.isTran=dis.isTran;
    
    if exist('FREE','var'), bcf.FREE=FREE; else bcf.FREE=0; end
    
    bcf.IBCFCB =MFLOWparvals(strmatchi('IBCFCB',MFLOWparnams),bcf.SCEN);
    bcf.HDRY   =MFLOWparvals(strmatchi('HDRY'  ,MFLOWparnams),bcf.SCEN);
    bcf.IWDFLG =MFLOWparvals(strmatchi('IWDFLG',MFLOWparnams),bcf.SCEN);
    bcf.WETFCT =MFLOWparvals(strmatchi('WETFCT',MFLOWparnams),bcf.SCEN);
    bcf.IWETIT =MFLOWparvals(strmatchi('IWETIT',MFLOWparnams),bcf.SCEN);
    bcf.IHDWET =MFLOWparvals(strmatchi('IHDWET',MFLOWparnams),bcf.SCEN);
    
    bcf.LAYCON =LAYparvals(:,strmatchi('LAYCON' ,LAYparnams));
    bcf.LAYAVG =LAYparvals(:,strmatchi('LAYAVG' ,LAYparnams));
    if GRID.AXIAL,
        if any(bcf.LAYCON>0)
            bcf.LAYAVG(:)=3; 
        else
            bcf.LAYAVG(:)=2;
        end
    end

    bcf.TPRY   =LAYparvals(:,strmatchi('TPRY'   ,LAYparnams));
    
    if exist('WETRY','var')
        bcf.WETRY=WETDRY;
    else
        bcf.WETDRY =LAYparvals(:,strmatchi('WETDRY' ,LAYparnams));
    end
    
    %% Check TRAN and HY using LAYCON
    % mfLab requires the following:
    %   if TRAN exists for all layers --> select according to LAYCON
    %   if HY   exists for all layers --> select according to LAYCON
    %   if size(TRAN,3)>nNonc & < nLay --> select the first nNonc
    %   if size(TRAN,3)>nConv & < nLay --> select the first nCoNv
    %   any other combination is an error.
    % This allows the user to always provide all TRAN layers (nLay) and as
    % many HY layers as to cover all variants he has in mind, without
    % altering anything but LAYCON.
    % 
    % In all cases, mf_setup below will select TRAN and HY to exactly match
    % the number of convertible and non-convertible layers in LAYCON.
    % Most user convenient is to specify TRAN for all layers and HY only
    % for the convertible layers always. Most of the time we need HY only
    % for the top most or the top couple of layers.
    % TO 120813
    
    % TRAN and HY can be dealt with separately, using LAYCON

    % Non-convertible layers (TRAN-layers)
    Inonc = find(bcf.LAYCON==0 | bcf.LAYCON==2);
    nNonc = numel(Inonc); % number of non convertible layers in LAYCON
    
    % Convertible layers (HY-layers)
    Iconv = find(bcf.LAYCON==1 | bcf.LAYCON==3);
    nConv = numel(Iconv);

    % TRAN first
    if exist('TRAN','var')
        if isvector(TRAN),
            TRAN = repmat(XS(TRAN(:)),[GRID.Ny,GRID.Nx,1]); % TRAN->3D
        end
        if size(TRAN,3) == GRID.Nlay,
            TRAN = TRAN(:,:,Inonc);    % select from full TRAN array
        elseif size(TRAN,3)>=nNonc && size(TRAN,3)<GRID.Nlay
            TRAN = TRAN(:,:,1:nNonc);  % select first nNonc layers
        else
            error('%s: Either specify TRAN (%d layers) for all %d layers or for the %d non convertable layers only, see LAYCON',...
                mfilename,size(TRAN,3),GRID.Nlay,nNonc);
        end
        if GRID.AXIAL
            bcf.TRAN = GRID.TWOPIR(:,:,1:size(TRAN,3)).*TRAN;
        else
            bcf.TRAN = TRAN;
        end
    end
    
    % Then HY
    if exist('HY'  ,'var')
        if isvector(HY),
            HY  =repmat(XS(HY(:))  ,[GRID.Ny,GRID.Nx,1]); % HY->3D
        end
        if size(HY,  3) == GRID.Nlay,
            HY = HY(:,:,Iconv);   % select from full HY array
        elseif size(HY,3)>=nConv && size(HY,3)<GRID.Nlay
            HY = HY(:,:,1:nConv); % select first nConv layers
        else
            error('%s: Either specify HY (%d layers) for all %d layers or for the %d non convertible layers only, see LAYCON',...
                mfilename,size(HY,3),GRID.Nlay,nConv);
        end
        if GRID.AXIAL
            bcf.HY = GRID.TWOPIR(:,:,1:size(HY,  3)).*HY;
        else
            bcf.HY = HY;
        end
    end        
        
    % VCONT
    if bcf.GRID.Nlay>1
        try            
            bcf.VCONT = VCONT;
            if GRID.AXIAL,
                bcf.VCONT = GRID.TWOPIR(:,:,1:size(bcf.VCONT,3)).*bcf.VCONT;
            end
        catch ME
            error('mfLab:mf_setup:VCONT_missing_in_workspace',...
                ['BCF pacakge is on, so VCONT must be specified in workspace if any LAYCBD>0 !\n',...
                'Maybe you must switch BCF on and LPF on in the NAM worksheet']);            
        end
    end
    if any(bcf.isTran)
        try
            bcf.SF1= SF1;
            if GRID.AXIAL,
                bcf.SF1 = GRID.TWOPIR(:,:,1:size(bcf.SF1,3)).*SF1;
            end            
        catch ME
            error('mfLab:mf_setup:SF1_missing_in_workspace',...
                'SF1 must be specified in workspace if any stress period is transient !');            
        end
        if any(bcf.LAYCON==2 | bcf.LAYCON==3)
            try
                bcf.SF2 = SF2;
                if GRID.AXIAL,
                    bcf.SF2 = GRID.TWOPIR(:,:,1:size(bcf.SF2,3)).*bcf.SF2;
                end
            catch ME
                error('mfLab:mf_setup:SF2_missing_in_workspace',...
                    'SF2 must be specified in workspace\nif any stress period is transient and any LAYCON=2 or 3 !');            
            end
        end
    end

    writeBCF(basename,bcf) 
end

%% ===== THE lpf-file ==================
lpf.SCEN=namSCEN('LPF',nam.PCKG,nam.SCEN);
upw.SCEN=namSCEN('UPW',nam.PCKG,nam.SCEN);
nwt     =namSCEN('NWT',nam.PCKG,nam.SCEN);

if lpf.SCEN || upw.SCEN

    lpf.FREE= bas.FREE;

    if lpf.SCEN
        lpf.name='LPF';
        lpf.upw = false;
        if upw.SCEN
            error(['%s: You can''t have LPF and UPW packages on together\n',...
                '   REMEDY: Check NAM sheet'],mfilename);
        end
        if nwt
            error(['%s: You can''t use packages LPF and NWT together\n',...
                   'REMEDY: Combine LPF with other solvers or combine UPW with NWT'],...
                   mfilename);
        end
    else
        lpf.name = 'UPW';        
        lpf.upw = true;
        if lpf.SCEN
            error(['%s: You can''t use LPF and UPW packages together\n',...
                   'REMEDY: Check NAM sheet'],mfilename);
        end
        if ~nwt
            error(['%s: You should use UPW and NWT together or not UPW and another solver.\n',...
                   'REMEDY: Combine UPW with NWT or LPF/BCF with other sovlers.'],mfilename);
        end
    end
        
    fprintf('Generating %s struct and file\n',lpf.name);
    lpf.unit=nam.UNIT(strmatchi(lpf.name,nam.PCKG));
    lpf.ext =nam.EXT {strmatchi(lpf.name,nam.PCKG)};
    lpf.GRID=GRID;

    lpf.I___CB =MFLOWparvals(strmatchi(['I' lpf.name 'CB'],MFLOWparnams),1); % file unit number for saving cell by cell flow terms   

    % Options for MF2005
    i=strmatchi('STORAGECOEFFICIENT'  ,MFLOWparnams,'exact');
    if i>0 && MFLOWparvals(i,1)~=0, lpf.STORAGECOEFFICIENT = MFLOWparvals(i,1); end

    i=strmatchi('CONSTANTCV'  ,MFLOWparnams,'exact');
    if i>0 && MFLOWparvals(i,1)~=0, lpf.CONSTANTCV = MFLOWparvals(i,1); end

    i=strmatchi('THICKSTRT'  ,MFLOWparnams,'exact');
    if i>0 && MFLOWparvals(i,1)~=0, lpf.THICKSTRT = MFLOWparvals(i,1); end

    i=strmatchi('NOCVCORRECTION'  ,MFLOWparnams,'exact');
    if i>0 && MFLOWparvals(i,1)~=0, lpf.NOCVCORRECTION = MFLOWparvals(i,1); end
    
    
    lpf.HDRY   =MFLOWparvals(strmatchi('HDRY',MFLOWparnams),1);
    lpf.NP___  =MFLOWparvals(strmatchi('NPLPF',MFLOWparnams),1);
    lpf.WETFCT =MFLOWparvals(strmatchi('WETFCT',MFLOWparnams),1);
    lpf.IWETIT =MFLOWparvals(strmatchi('IWETIT',MFLOWparnams),1);
    lpf.IHDWET =MFLOWparvals(strmatchi('IHDWET',MFLOWparnams),1);

    if lpf.upw
    % IPHDRY is a flag that indicates whether groundwater head will be set to HDRY
    % when the groundwater head is less than 1?10-4 above the cell bottom (units
    % defined by LENUNI). If IPHDRY=0, then head will not be set to HDRY.
    % If IPHDRY>0, then head will be set to HDRY.
        lpf.IPHDRY = MFLOWparvals(strmatchi('IPHDRY',MFLOWparnams),1);
    end
    
    lpf.LAYCON =LAYparvals(:,strmatchi('LAYCON' ,LAYparnams)); % Layer type 1=convertible
        
    lpf.LAYAVG =LAYparvals(:,strmatchi('LAYAVG' ,LAYparnams)); % Hor cond comp method (0=harmonic)
    if GRID.AXIAL
        if any(lpf.LAYCON)>0
            lpf.LAYAVG(:)=2;  % most appropriate for unconfined flow
        else 
            lpf.LAYAVG(:)=1;  % most appropriate for confined flow
        end
    end
    
    lpf.LAYWET =LAYparvals(:,strmatchi('LAYWET' ,LAYparnams)); % Layer wettability flag
    if lpf.upw  % LAYWET must be zeros when upw is used
        lpf.LAYWET(:)=0;
    end
    
    if any(GRID.LAYCBD>0)
        try
            lpf.VKCB   =VKCB; % vertical conductivity of confining units must be specified in matlab
            if GRID.AXIAL
                lpf.VKCB=GRID.TWOPIR(:,:,1:gr.Ncbd).*lpf.VKCB;
            end
        catch ME
            error('mfLab:mf_setup:VKCB_missing_in_workspace',...
                'VKCB must be defined in mf_adapt if any LAYCBD~=0 !');
        end
    end
    
    lpf.WETDRY =LAYparvals(:,strmatchi('WETDRY' ,LAYparnams)); % Layer wetting method flag
    
    % allow specifiation of WETDRY as cell array also in mf_adapt
    if exist('WETDRY','var'),     % if not, we keep bcf.WETDRY as is
        if ~iscell(WETDRY),
            error('WETDRY specified in mf_adapt must be a cell array with one model layer WETDRY values per array cell');
        end
        wetdry=lpf.WETDRY;        % store WETDRY from spreadsheet LAY
        lpf.WETDRY=cell(size(WETDRY));
        for iL=1:size(lpf.WETDRY,1)  % may be smaller than number of layers
            if wetdry(iL)==0,        % this triggers use of mf_adapt values
                lpf.WETDRY{iL}=WETDRY{iL};  % default
            else
                lpf.WETDRY{iL}=wetdry(iL);
            end
        end
    end            

    % LAYVKA in workspace takes precedence over LAYVKA in workbook/lay
    if exist('LAYVKA','var'),        
        lpf.LAYVKA = LAYVKA;
    else
        lpf.LAYVKA =XS(LAYparvals(:,strmatchi('LAYVKA' ,LAYparnams)));
    end
        
    % VK (always defined in workspace) takes precedence over LAYVKA
    if exist('VK','var')        
        lpf.LAYVKA(:)=0; % meaning VKA is VK
        lpf.VKA = VK;
        if GRID.AXIAL
            lpf.VKA=GRID.TWOPIR.*lpf.VKA;
        end
    else     % VKA must be specified in the workspace
        try
            lpf.VKA = VKA;  % VKA could be given as a column vector
            if ismatrix(VKA) && GRID.Nlay>1
                lpf.VKA=XS(VKA(:));
            end
            if GRID.AXIAL
                lpf.VKA=GRID.TWOPIR.*lpf.VKA;
            end
        catch ME
            error('mfLab:mf_setup:VKA_missing_in_workspace',...
                'You must define VKA or VK in mf_adapt if you use %s package!',lpf.name);
        end
    end
    
    % horizontal conductivity in x-direction
    try
        lpf.HK  = HK;
        if GRID.AXIAL,
            lpf.HK =GRID.TWOPIR.*lpf.HK;
        end
    catch ME
        error('mfLab:mf_setup:HK_missing_in_workspace',...
            'You must define HK in mf_adapt if you use the %s package!',lpf.name);
    end
    
    % CHANI in workspace takes precedence over CHANI in workbook/LAY
    if exist('CHANI','var')
        lpf.CHANI=XS(CHANI(:));
    else
        lpf.CHANI  =XS(LAYparvals(:,strmatchi('CHANI'  ,LAYparnams))); % Layer wetting method flag
    end
    
    % HANI must be specified in the workspace if any lpf.CHANI<=0
    if any(lpf.CHANI<=0)
        try
            lpf.HANI=HANI; % HANI could be given as full 3D array
            if ismatrix(HANI) && GRID.Nlay>1, % if not make it 3D
                lpf.HANI=XS(HANI(:));
            end
        catch ME
            error('mf_setup:lpfupw:HANI_missing_in_workspace',...
            'if any CHANI<=0, HANI must be specified in mf_adapt!');
        end
    end
    
    lpf.isTran = dis.isTran;  % whether or not transient flow
    
    if any(lpf.isTran)
        try
            lpf.SS = SS;
            if GRID.AXIAL,
                lpf.SS=GRID.TWOPIR(:,:,1:size(lpf.SS,3)).*lpf.SS;
            end
        catch ME
            error('mf_setup:lpwupw:SS_missing_in_workspace',...
            'if transient SS must be specified in mf_adapt!');
        end
        if any(lpf.LAYCON),
            try
                lpf.SY = SY;
                if GRID.AXIAL,
                    lpf.SY=GRID.TWOPIR(:,:,1:size(lpf.SY,3)).*SY;
                end
            catch ME
                error('mf_setup:lpfupw:SY_missing_in_workspace',...
                'if transient SY must be specified in mf_adapt!');        
            end
        end
    end
   
    lpf.LAYWET(lpf.LAYCON==0)=0; % LAYWET must be 0 if LAYCON is 0

    writeLPF(basename,lpf); %writs both lpf and upw
end

%% ===== the VDF-file for SEAWAT ================================
vdf.SCEN=namSCEN('VDF',nam.PCKG,nam.SCEN);

if vdf.SCEN>0
    fprintf('Generating VDF struct and file\n');
    vdf.unit=nam.UNIT(strmatchi('VDF',nam.PCKG));
    vdf.ext =nam.EXT {strmatchi('VDF',nam.PCKG)};
    vdf.FREE=bas.FREE;
    vdf.NLAY=GRID.Nlay;
    vdf.NPER=NPER;
    
    % for each simulation
    %1  Note alternative parameter names depending on version
    if strmatchi('MT3DRHOFLG',SEAWATparnams)
        vdf.MT3DRHOFLG=SEAWATparvals(strmatchi('MT3DRHOFLG',SEAWATparnams),vdf.SCEN); % recharge option, in what layer
    else
        vdf.MT3DRHOFLG=SEAWATparvals(strmatchi('MTDNCONC'   ,SEAWATparnams),vdf.SCEN); % recharge option, in what layer
    end
    vdf.MFNADVFD=SEAWATparvals(strmatchi('MFNADVFD',SEAWATparnams),vdf.SCEN); 
    vdf.NSWTCPL =SEAWATparvals(strmatchi('NSWTCPL' ,SEAWATparnams),vdf.SCEN); 
    vdf.IWTABLE =SEAWATparvals(strmatchi('IWTABLE' ,SEAWATparnams),vdf.SCEN); 

    %2
    vdf.DENSEMIN =SEAWATparvals(strmatchi('DENSEMIN',SEAWATparnams),vdf.SCEN); 
    vdf.DENSEMAX =SEAWATparvals(strmatchi('DENSEMAX',SEAWATparnams),vdf.SCEN); 
 
    %3
    vdf.DNSCRIT =SEAWATparvals(strmatchi('DNSCRIT' ,SEAWATparnams),vdf.SCEN);
    
    if vdf.MT3DRHOFLG>=0
        %4
        vdf.DENSEREF=SEAWATparvals(strmatchi('DENSEREF',SEAWATparnams),vdf.SCEN); 
        vdf.DRHODC  =SEAWATparvals(strmatchi('DRHODC(1)',SEAWATparnams),vdf.SCEN);
        if vdf.MT3DRHOFLG==0
            %6
            vdf.INDENSE =PERvals(:,strmatchi('INDENSE',PERnams)); % if >0 read dense for each period
            if any(vdf.INDENSE)>0
                %7
                vdf.DENSE=DENSE;  % get DENS{NPER}(NROW,NCOL) from mf_adapt
            end
        end

    elseif vdf.MT3DRHOFLG==-1
    %4a
        vdf.DENSEREF  =SEAWATparvals(strmatchi('DENSEREF' ,SEAWATparnams),vdf.SCEN);
        vdf.DRHODPRHD =SEAWATparvals(strmatchi('DRHODPRHD',SEAWATparnams),vdf.SCEN); 
        vdf.PRHDREF   =SEAWATparvals(strmatchi('PRHDREF'  ,SEAWATparnams),vdf.SCEN); 
        %4b
        vdf.NSRHOEOS  =SEAWATparvals(strmatchi('NSRHOEOS' ,SEAWATparnams),vdf.SCEN);
        
        NCOMP=MT3Dvals(strmatchi('NCOMP',MT3Dnams),1);
        if vdf.NSRHOEOS>NCOMP
            warning('on','mf_setup:VDF:NrOfSpeciesVDFnotequaltoNCOMP');
            warning('mf_setup:VDF:NrOfSpeciesVDFnotequaltoNCOMP',...
                ['The number of species NSRhoEOS=%d in VDF that are used for the density computation\n',...
                 'is larger than the number of species defined in BTN (NCOMP=%d).\n',...
                'mfLab sets NSRHOEOS==NCOMP.\n'],vdf.NSRHOEOS,NCOMP);        
            vdf.NSRHOEOS = NCOMP;
        end
        
        %4c
        vdf.MTRHOSPEC  =SEAWATparvals(strmatchi('MTRHOSPEC',SEAWATparnams),1:vdf.NSRHOEOS);
        vdf.DRHODC     =SEAWATparvals(strmatchi('DRHODC'   ,SEAWATparnams,'exact'),1:vdf.NSRHOEOS); 
        vdf.CRHOREF    =SEAWATparvals(strmatchi('CRHOREF'  ,SEAWATparnams),1:vdf.NSRHOEOS);        
    end
    %5
    vdf.FIRSTDT =SEAWATparvals(strmatchi('FIRSTDT' ,SEAWATparnams),vdf.SCEN); 
     
    writeVDF(basename,vdf) 
end

%% ===== the VSC-file for SEAWAT ================================
vsc.SCEN=namSCEN('VSC',nam.PCKG,nam.SCEN);

if vsc.SCEN
    fprintf('Generating VSC struct and file\n');
    vsc.unit=nam.UNIT(strmatchi('VSC',nam.PCKG));
    vsc.ext =nam.EXT {strmatchi('VSC',nam.PCKG)};
    vsc.FREE=bas.FREE;
    vsc.NLAY=GRID.Nlay;
    vsc.NPER=NPER;
    
    if ~strmatchi('VDF',nam.PCKG),
        error('mfLab:mf_setup:vscPackageRequiresVDFpackage',...
            ['VSC package only works if VDF package is also on, see swt_v4 manual page 20,\n',...
             'section: ''Iput instructions and Evaluation of Temperature Output''']);
    end
    
    % for each simulation
    %1
    vsc.MT3DMUFLG=SEAWATparvals(strmatchi('MT3DMUFLG',SEAWATparnams),vsc.SCEN); % recharge option, in what layer

    %2
    vsc.VISCMIN  =SEAWATparvals(strmatchi('VISCMIN'  ,SEAWATparnams),vsc.SCEN); 
    vsc.VISCMAX  =SEAWATparvals(strmatchi('VISCMAX'  ,SEAWATparnams),vsc.SCEN); 
    vsc.VISCREF  =SEAWATparvals(strmatchi('VISCREF'  ,SEAWATparnams),vsc.SCEN); 
    
    vsc.NSMUEOS  =SEAWATparvals(strmatchi('NSMUEOS'  ,SEAWATparnams),vsc.SCEN); 
    %3 
    vsc.MTMUSPEC =SEAWATparvals(strmatchi('MTMUSPEC' ,SEAWATparnams),:);
    vsc.DMUDC    =SEAWATparvals(strmatchi('DMUDC'    ,SEAWATparnams),:);
    vsc.CMUREF   =SEAWATparvals(strmatchi('CMUREF'   ,SEAWATparnams),:);
    
    vsc.MTMUTEMPSPEC=SEAWATparvals(strmatchi('MTMUTEMPSPEC',SEAWATparnams),vsc.SCEN);
    vsc.MUTEMPOPT   =SEAWATparvals(strmatchi('MUTEMPOPT',SEAWATparnams),vsc.SCEN);  
    vsc.AMUCOEF     =SEAWATparvals(strmatchi('AMUCOEF'     ,SEAWATparnams),vsc.MUTEMPOPT);
    vsc.AMUCOEF(isnan(vsc.AMUCOEF))=[];

        
    if vsc.MT3DMUFLG==0
        %4
        vsc.INVISC = PERvals(:,strmatchi('INVISC',   PERnams));
        if any(vsc.INVISC)>0
            %5
            vsc.VISC=VISC; % get VISC{NPER}(NROW,NCOL) from mf_adapt
        end
    end


    writeVSC(basename,vsc) 
end

%% ===== THE rch-file ============================
rch.SCEN=namSCEN('RCH',nam.PCKG,nam.SCEN);

if rch.SCEN
    fprintf('Generating RCH struct and file\n');
    rch.unit=nam.UNIT(strmatchi('RCH',nam.PCKG));
    rch.ext =nam.EXT {strmatchi('RCH',nam.PCKG)};
    
    rch.FREE=bas.FREE;
    rch.NPER=NPER;
    % for each simulation
    %1
    rch.NPRCH=0;         % we will use no parameters
    %2    
    rch.NRCHOP =MFLOWparvals(strmatchi('NRCHOP',MFLOWparnams),rch.SCEN); % recharge option, in what layer
    rch.IRCHCB =MFLOWparvals(strmatchi('IRCHCB',MFLOWparnams),rch.SCEN); % file unit number for saving cell by cell flow terms
    
    %3
    %4
    
    %% for each stress period
    %5  
    
    rch.INRECH =PERvals(:,strmatchi('INRECH',PERnams)); % if >0 read recharge layer
    rch.INIRCH =PERvals(:,strmatchi('INIRCH',PERnams)); % skipped unless recharge in specific layers

    rch.rech   =PERvals(:,strmatchi('RECH'  ,PERnams));
    rch.irch   =PERvals(:,strmatchi('IRCH'  ,PERnams));
    
    %6
    if any(rch.INRECH==0) && ~exist('RECH','var'); % RECH must be defined in mf_adapt
        error('If any IN%s in the PER worksheet = 0 then you must define %s in mf_adapt !','RECH','RECH');
    end
    
    % if RECH exists as a list, make it 3D with 1:NPER into 3rd dimension
    if exist('RECH','var') && ~iscell(RECH) && size(RECH,3)<NPER
        nper=size(RECH,3);
        RECH(:,:,NPER)=NaN; % Allocate by extending RECH
        if all(size(RECH(:,:,1))==1) || all(size(RECH(:,:,1))==gr.size(1:2))
            % extend RECH to a 3D array [Ny,Nx,NPER]
            RECH(:,:,size(RECH,3)+1:NPER)= ...
                bsxfun(@times,RECH(:,:,end),ones(1,1,NPER-size(RECH,3)));
        else
            error('mf_setup:RECH:WrongDimensions',...
                ['%s; first 2 dimensions of RECH =[%d,%d] must match [Ny,Nx]=[%d,%d] or be 1,1\n',...
                 'Make sure you provide RECH as 3D (Ny,Nx,NPER) array or (1,1,NPER) array,\n',...
                 'to allow mfLab to recognize time from space in the RECH array.'],...
                size(RECH(:,:,1)),gr.size(1:2));
        end            
    end
        
    % Handle axial symmetry
    for iPer=1:rch.NPER
        if rch.INRECH(iPer)==0
            if iscell(RECH)
                if GRID.AXIAL
                    rch.RECH{iPer}=GRID.TWOPIR(:,:,1).*RECH{iPer};
                else
                    rch.RECH{iPer}=RECH{iPer};
                end
            else
                if GRID.AXIAL 
                    rch.RECH{iPer}=GRID.TWOPIR(:,:,1).*RECH(:,:,iPer);
                else
                    rch.RECH{iPer}=RECH(:,:,iPer);
                end
            end
        else
            if GRID.AXIAL 
                rch.RECH{iPer}=GRID.TWOPIR(:,:,1).*rch.rech(iPer);
            else
                rch.RECH{iPer}=rch.rech(iPer);
            end
        end
    end

    %7
    %8
    if any(rch.INIRCH==0) && ~exist('IRCH','var')
        error('If any IN%s in the PER worksheet = 0 then you must define %s in mf_adapt !',...
                'IRCH','IRCH');
    end
    for iPer=1:rch.NPER
        if rch.INIRCH(iPer)==0
            if iscell(IRCH)
                rch.IRCH{iPer}=IRCH{iPer};
            else
                %% Allow IRCH to given for only the first stress period
                rch.IRCH{iPer}=IRCH(:,:,min(size(IRCH,3),iPer));
            end
        else
            rch.IRCH{iPer}=rch.irch(iPer);
        end
    end
    
    writeRCH(basename,rch) 
end

%% s===== THE evt-file (the Evaporation package) ==============
evt.SCEN=namSCEN('EVT',nam.PCKG,nam.SCEN);
ets_SCEN=namSCEN('ETS',nam.PCKG,nam.SCEN); % Note the underscore here

% we use struct evt for both the packages evt and ets because the code has
% 95% overlap and both packages are mutually exclusive
% if ETS is on, columns ISGDF and PXDP.. and PETM.. are read, where .. is any
% sequence, so that multiple columns are used to defined the segments of
% the extinction curve.

if evt.SCEN || ets_SCEN
    if evt.SCEN,
        fprintf('Generating EVT struct and file\n');
        evt.unit=nam.UNIT(strmatchi('EVT',nam.PCKG));
        evt.ext =nam.EXT {strmatchi('EVT',nam.PCKG)};
    end
    if ets_SCEN,
        fprintf('Generating ETS struct and file\n');
        evt.unit=nam.UNIT(strmatchi('ETS',nam.PCKG));
        evt.ext =nam.EXT {strmatchi('ETS',nam.PCKG)};
    end
    
    evt.FREE=bas.FREE;
    evt.NPER=NPER;
    
    evtp.NPEVT=0;
    %1
    %2
    evt.NEVTOP =MFLOWparvals(strmatchi('NEVTOP',MFLOWparnams),1);
    evt.NEVTOP =max(min(evt.NEVTOP,2),1);
    evt.IEVTCB =MFLOWparvals(strmatchi('IEVTCB',MFLOWparnams),1); % file unit number for saving cell by cell flow terms
    %3
    %4
    %5
    evt.INSURF =PERvals(:,strmatchi('INSURF',PERnams));
    evt.INEVTR =PERvals(:,strmatchi('INEVTR',PERnams));
    evt.INEXDP =PERvals(:,strmatchi('INEXDP',PERnams));
    evt.INIEVT =PERvals(:,strmatchi('INIEVT',PERnams));
    
    evt.surf   =PERvals(:,strmatchi('SURF',PERnams));
    evt.evtr   =PERvals(:,strmatchi('EVTR'  ,PERnams));
    evt.exdp   =PERvals(:,strmatchi('EXDP'  ,PERnams));
    evt.ievt   =PERvals(:,strmatchi('IEVT'  ,PERnams));
    
    if ets_SCEN
        evt.pxdp   =PERvals(:,strmatchi('PXDP'  ,PERnams));
        evt.petm   =PERvals(:,strmatchi('PETM'  ,PERnams));
    end
    evt.SURF   =cell(evt.NPER,1); for i=1:evt.NPER, evt.SURF{i}=evt.surf(i); end
    evt.EVTR   =cell(evt.NPER,1); for i=1:evt.NPER, evt.EVTR{i}=evt.evtr(i); end
    evt.EXDP   =cell(evt.NPER,1); for i=1:evt.NPER, evt.EXDP{i}=evt.exdp(i); end
    evt.IEVT   =cell(evt.NPER,1); for i=1:evt.NPER, evt.IEVT{i}=evt.ievt(i); end

    if ets_SCEN  % extra for ets, read the columns PXDP and PETM

        evt.INSGDF =PERvals(:,strmatchi('INSGDF',PERnams)); % segment definition read flag

        IPXDP      =strmatchi('PXDP'  ,PERnams,'exact');
        IPETM      =strmatchi('PETM'  ,PERnams,'exact');
        
        if any(IPXDP)==0 || any(IPETM)==0
            evt.NETSEG=1;       % in this case package ETS == package EVT
        else
        
            if numel(IPETM) ~= numel(IPXDP),
                error('Number of PXDP columns (%d) must equal number of PETM columns (%d) in worksheet PER',...
                    numel(IPXDP),numel(IPETM));
            end
        
            % The umber of sections is 1+the number of PXDP or PETM columns in
            % the PER worksheet ! NETSEG==1 implies ETS package == EVT package
            evt.NETSEG=numel(IPXDP)+1;
    
            evt.PXDP   =cell(evt.NPER,evt.NETSEG-1);
            evt.PETM   =cell(evt.NPER,evt.NETSEG-1);
        end
    end
       
    %6
    % if any INSURF in the spreadsheet==0, the SURF must exist in Matlab's
    % workspace at this point
    if any(evt.INSURF==0) && ~exist('SURF','var')
                   error('If any IN%s in the PER worksheet = 0 then you must define %s in mf_adapt !','SURF','SURF');
    end
    
    % SURF must always be defined for the first stress period
    % See comment in worksheet
    % if INSURF(1)==-2, SURF(1) is interpreted as m above ground surface
    %   so use negative value in sheet defining distance of SURF below
    %   ground surface
    if evt.INSURF(1)==-2; evt.SURF{1}=evt.surf(1)+GRID.Z(:,:,1); evt.INSURF(1)=2; end
    % if INSURF(1)==-1, then SURF(1) is interprete as uniform SURF elevation
    if evt.INSURF(1)==-1; evt.SURF{1}=evt.surf(1);               evt.INSURF(1)=1; end
    
    for iPer=1:evt.NPER
        if evt.INSURF(iPer)==0
            if iscell(SURF)
                evt.SURF{iPer}=SURF{iPer};
            else
                % use INSURF(iPER)<0 for IPER>0 to repeat previous values
                % then SURF(iPER) does not matter.
                evt.SURF{iPer}=SURF(:,:,iPer);
            end
        end
    end
    
    %7
    if any(evt.INEVTR==0) && ~exist('EVTR','var')
        error('If any IN%s in the PER worksheet = 0 then you must define %s in mf_adapt !','EVTR','EVTR');
    end
    
     % if EVTR exists as a list, make it 3D with 1:NPER into 3rd dimension
    if exist('EVTR','var') && ~iscell(EVTR) && any(size(EVTR)==numel(EVTR)),
        EVTR=permute(EVTR(:),[2,3,1]);
    end

    for iPer=1:evt.NPER
        if evt.INEVTR(iPer)==0
            if iscell(EVTR)
                if GRID.AXIAL, evt.EVTR{iPer} =GRID.TWOPIR(:,:,1).*EVTR{iPer}; else evt.EVTR{iPer} =EVTR{iPer}; end
            else
                if GRID.AXIAL, evt.EVTR{iPer} =GRID.TWOPIR(:,:,1).*EVTR(:,:,iPer); else evt.EVTR{iPer} =EVTR(:,:,iPer); end
            end
        else
            if GRID.AXIAL, evt.EVTR{iPer}=GRID.TWOPIR(:,:,1).*evt.evtr(iPer); else evt.EVTR{iPer}=evt.evtr(iPer); end
        end
    end
    
    %8 parameters not yet implemented
    
    %9
    if any(evt.INEXDP==0) && ~exist('EXDP','var'); % EVTR must by defined in mf_adapt
        error('If any IN%s in the PER worksheet = 0 then you must define %s in mf_adapt !','EXDP','EXDP');
    end
    if evt.INEXDP(1)==-1; evt.EXDP{1}=evt.exdp(1); evt.INEXDP(1)=1; end
    for iPer=1:evt.NPER
        if evt.INEXDP(iPer)==0
            if iscell(EXDP)
                evt.EXDP{iPer}=EXDP{iPer};
            else
                evt.EXDP{iPer}=EXDP(:,:,iPer);
            end
        end
    end
    
    %10
    if evt.NEVTOP>1
        if any(evt.INIEVT==0) && ~exist('IEVT','var'); % EVTR must by defined in mf_adapt
            error('If any IN%s in the PER worksheet = 0 then you must define %s in mf_adapt !','IEVT','IEVT');
        end
        for iPer=1:evt.NPER
            if evt.INIEVT(iPer)==0
                if iscell(IEVT)
                    evt.IEVT{iPer}=IEVT(iPer);
                else
                    evt.IEVT{iPer}=IEVT(:,:,iPer);
                end
            end
        end
    end    

    
    if evt.SCEN % ETS package, watch underscore, ets_SCEN is a scalar boolean
        writeEVT(basename,evt,evtp) 
    end
    
    if ets_SCEN % ETS package        
        if evt.NETSEG>1
            if any(evt.INSGDF==0) && ~exist('PXDP','var'); % PXDP must by defined in mf_adapt
                error('If any ISGDF in the PER worksheet == 0 then you must define %s in mf_adapt !','PXDP');
            end
            if any(evt.INSGDF==0) && ~exist('PETM','var') % PETM must be defined in mf_adapt
                error('If any ISGDF in the PER worksheet == 0, then you must define %s in mf_adapt !','PETM');
            end
        
            for iPer=1:evt.NPER
                 if evt.INSGDF(iPer)==0
                     if iscell(PXDP)
                         evt.PXDP{iPer,:}=PXDP(iPer,:);
                     else
                         evt.PXDP{iPer}=PXDP(:,:,iPer);
                     end
                     if iscell(PETM)
                        evt.PETM{iPer,:}=PETM(iPer,:);
                     else
                         evt.PETM{iPer}=PETM(:,:,iPer);
                     end
                 else
                     evt.PXDP{iPer}=evt.pxdp(iPer);
                     evt.PETM{iPer}=evt.petm(iPer);
                 end
            end
        end
        
        writeETS(basename,evt,evtp) 
    end

end

%% ===== MNW1 multi-node well package (2002) ====================================
mnw1.SCEN=namSCEN('MNW1',nam.PCKG,nam.SCEN);  % check if this boundary package is on

% this code is only used to generate the MNW file directly from mnw1
% objects without first generating the MNW and PNTSRC lists. Hence it is
% expected that MNW1 exists, the MNW1 package is on and the class of MNW1 is
% MNW1Obj.

if mnw1.SCEN  % if this boundary package is "on" then ...
    fprintf('Generating %s struct and file\n','MNW1');
        
    % UNIT numbers for cell-by-cell flow term output files, specified in
    
    mnw1.FREE = bas.FREE;
    
    % negative mnw1.ICB makes sure intermedate data is printed in list file
    mnw1.ICB      = -MFLOWparvals(strmatchi('IWELCB', MFLOWparnams),mnw1.SCEN);
    
    mnw1.IWELPT   = MFLOWparvals(strmatchi('IWELPT', MFLOWparnams),mnw1.SCEN);
    mnw1.WELL1    = MFLOWparvals(strmatchi('WELL1' , MFLOWparnams),mnw1.SCEN);
    mnw1.BYNODE   = MFLOWparvals(strmatchi('BYNODE', MFLOWparnams),mnw1.SCEN);
    mnw1.QSUM     = MFLOWparvals(strmatchi('QSUM'  , MFLOWparnams),mnw1.SCEN);
    mnw1.LOSSTYPE = MFLOWparvals(strmatchi('LOSSTYPE', MFLOWparnams),mnw1.SCEN);
  
    idx=strmatchi('MNW1',nam.PCKG);
    mnw1.UNIT  = nam.UNIT(idx); % unit number of input file of this package
    mnw1.allUNITS = nam.UNIT; % all unit numbers in use
    mnw1.ext =nam.EXT {idx}; % extension of input file of this package

    mnw1.NPER =NPER;

    %% Collect all MNW1Obj into mnw1.MNW
    I=strmatchi('MNW1Obj',{mf_variables.class},'exact');
    if I
        if ~isfield(mnw1,'MNW'), mnw1.MNW=[]; end
        for i=1:numel(I)
            eval(['mnw1.MNW = [mnw1.MNW ; ',mf_variables(I(i)).name '];']);
            warning('on','mf_setup:BCNtypeWell:MNW1ObjAddedAutomatically');
            warning('mf_setup:BCNtypeWell:MNW1ObjAddedAutomatically',...
                '%s: Generating BCNtype = %s: multinode wells (type 1) <<%s>> of class <<%s>> added',...
                mfilename,'MNW',mf_variables(I(1)).name,mf_variables(I(1)).class);
        end
    else
        error(['%s: MNW1 is switched on in the NAM sheet, but not MNW1Obj objects found in workspace.\n',...
            'REMEDY: you might have MNW2Obj objects or wellObj objects instead.\n',...
            '        in that case switch on the correct packages in the NAM sheet'],...
            mfilename);
    end
    
    MXSS = writeMNW1(basename,mnw1);
end
clear('mnw1'); % and do away with mnw1

%% ===== MNW2 multinode pakcage (2009) ====================================
mnw2.SCEN=namSCEN('MNW2',nam.PCKG,nam.SCEN);  % check if this boundary package is on

if mnw2.SCEN  % if this boundary package is "on" then ...
    fprintf('Generating %s struct and file\n','MNW2');

    % UNIT numbers for cell-by-cell flow term output files, specified in
    mnw2.ICB      = MFLOWparvals(strmatchi('IWELCB',   MFLOWparnams),mnw2.SCEN);
    mnw2.MNWPRNT  = MFLOWparvals(strmatchi('MNWPRNT',  MFLOWparnams),mnw2.SCEN);
    mnw2.MNW      = MNW;  % array of mnwobj 
  

    idx=strmatchi('MNW2',nam.PCKG);
    mnw2.unit=nam.UNIT(idx); % unit number of input file of this package
    mnw2.ext =nam.EXT {idx}; % extention of input file of this package

    mnw2.FREE = bas.FREE;

    mnw2.NPER = NPER;

    %% Collect all MNW2Obj into mnw2.MNW
    I=strmatchi('MNW2Obj',{mf_variables.class},'exact');
    if I
        if ~isfield(mnw2,'MNW'), mnw2.MNW=[]; end
        for i=1:numel(I)
            eval(['mnw2.MNW = [mnw2.MNW ; ',mf_variables(I(i)).name '];']);
            warning('on','mf_setup:BCNtypeWell:MNW2ObjAddedAutomatically');
            warning('mf_setup:BCNtypeWell:MNW2ObjAddedAutomatically',...
                '%s: Generating BCNtype = %s: multinode wells (type 1) <<%s>> of class <<%s>> added',...
                mfilename,'MNW',mf_variables(I(1)).name,mf_variables(I(1)).class);
        end
    else
        error(['%s: MNW2 is switched on in the NAM sheet, but not MNW2Obj objects found in workspace.\n',...
            'REMEDY: you might have MNW2Obj objects or wellObj objects instead.\n',...
            '        in that case switch on the correct packages in the NAM sheet'],...
            mfilename);
    end

    writeMNW2(basename,mnw2);
end
clear('mnw2'); % and do away with mnw2

%% ===== THE STRESSES WEL, GHB, DRN, RIV, CHD ===================================

IO = strmatchi('pointObj',{mf_variables.class},'empty');
pointObjects = eval(['[' sprintf(' %s',mf_variables(IO).name) ' ];']);
IO = strmatchi('lineObj',{mf_variables.class},'empty');
lineObjects = eval(['[' sprintf(' %s',mf_variables(IO).name) ' ];']);
IO = strmatchi('area2Obj',{mf_variables.class},'empty');
areaObjects = eval(['[' sprintf(' %s',mf_variables(IO).name) ' ];']);

% TODO extend with new types such as DRT  TO120903
bctypes={'WEL','GHB','DRN','RIV','CHD'};

activeStresses = bctypes(ismember(bctypes,nam.PCKG));

for ibc=1:numel(activeStresses)  % for each BCN type one loop
        
    bcn.FREE = bas.FREE;
    bcn.gr   = GRID;
    
    bcn.type   = activeStresses{ibc};
        
    idx=strmatchi(bcn.type,nam.PCKG,'exact');
    bcn.unit=nam.UNIT(idx); % unit number of input file of this package
    bcn.ext =nam.EXT {idx}; % extention of input file of this package
    
    I = strmatchi(bcn.type,MFLOWtxtnams,'once');
    if I
        bcn.AUX = MFLOWtxtvals(I,:);
        bcn.AUX = bcn.AUX(cellfun(@ischar,bcn.AUX));
    end
    
    if ~strmatchi(bcn.type,who,'exact')
        msgId = 'setup:NoDataInWorkspace';
        warning('on',msgId);
        warning(msgId,'No array <<%s>> found in workspace',bcn.type);
    else
        bcn.BCN = eval(bcn.type);
    end
        
    fprintf('%s: Generating %s struct and file\n',mfilename,bcn.type);
    
    % add line and area objects for the current stress type
    % we need separate bcType to deal with FLUX and WEL equivalents that
    % are allowed in pointObj, lineObj and area2Obj
    if strcmpi(bcn.type,'WEL')
        bcType={'FLUX','WEL'};
    else
        bcType=bcn.type;
    end
       
    if exist('pointObjects','var') && ~isempty(pointObjects)
        I  = strmatchi(bcType,{pointObjects.type},'empty');
        if ~isempty(I)
            bcn.point = pointObjects(I);
            bcn.point = bcn.point.check(IBOUND);
            % pointObjects(I)=[];  % SSM needs them
        end
    end
    if exist('lineObjects','var') && ~isempty(lineObjects)
        I  = strmatchi(bcType,{lineObjects.type},'empty');
        if ~isempty(I)
            bcn.line = lineObjects(I);
            bcn.line = bcn.line.check(IBOUND);
            % lineObjects(I)=[]; % SSM needs them
        end
    end
    if exist('areaObjects','var') && ~isempty(areaObjects);
        I  = strmatchi(bcType,{areaObjects.type},'empty');
        if ~isempty(I),
            bcn.area = areaObjects(I);
            bcn.area = bcn.area.check(IBOUND);
            % areaObjects(I)=[]; % SSM needs them
        end
    end
    % at this point we don't need bcType anymore and continue with bcn.type
    % which has WEL instead of FLUX which is allowed in pointObj, lineObj
    % and areaObj.

    % UNIT numbers for cell-by-cell flow term output files, specified in
    % worksheet MFLOW
    switch bcn.type
        case 'WEL'
            % unit number for basename.WEL file
            bcn.ICB =MFLOWparvals(strmatchi('IWELCB',MFLOWparnams),1);

            % get wellObj's from workspace
            I=strmatchi('wellObj',{mf_variables.class});

            if I  % workspace has wellObj

                % if first round, then add field to bcn struct
                if ~isfield(bcn,'well'), bcn.well=[]; end

                % for all wellObj's in the workspace ...
                for i=1:numel(I)

                    % add it to existing list in bcn struct
                    eval(['bcn.well = [bcn.well ; ',mf_variables(I(i)).name '];']);
                    fprintf('%s: Generating BCNtype = %s: wells <<%s>> of class <<%s>> added',...
                        mfilename,bcn.type,mf_variables(I(1)).name,mf_variables(I(1)).class);
                end
            end

        case 'GHB'
            bcn.ICB =MFLOWparvals(strmatchi('IGHBCB',MFLOWparnams),1);
        case 'DRN',
            bcn.ICB =MFLOWparvals(strmatchi('IDRNCB',MFLOWparnams),1);                
        case 'RIV',
            bcn.ICB =MFLOWparvals(strmatchi('IRIVCB',MFLOWparnams),1);            
        case 'CHD',
            bcn.ICB =MFLOWparvals(strmatchi('ICHDCB',MFLOWparnams),1);

            % verify number of CHD cols in case VDF is on because of
            % auxilliary arguments TO 130820
            if ismember('VDF',nam.PCKG)
                                
                % possible seawat options
                if exist('CHDDENSOPT','var')
                    bcn.AUX{1}='CHDDENSOPT';
                else
                    warning('on','mf_setup:CHD:CHDDENSTOPT');
                    warning('mf_setup:CHD:CHDDENSTOPT','CHDDENSTOP is active.\n');
                end
                if exist('CHDDEN'   ,'var')
                    % if BCN=CHD has 8 columns CHDDEN auxiiary variable must be set
                    % We check and do this here as service to the innocent user
                    % Only changes if CHDDEN opt is changed by USGS
                    bcn.AUX{2}='CHDDEN'   ;
                end

                if ~isfield(bcn,'AUX'),
                      nAux=0;
                else  nAux=sum(~cellfun('isempty',bcn.AUX));
                end
                
                %if only pointObj, lineObj and area2Obj are used, BCN does
                %not have to exist in the workspace
                if isfield(bcn,'BCN')
                    if iscell(bcn.BCN), ncol=size(bcn.BCN{1},2); else ncol=size(bcn.BCN,2); end

                    if ncol<6, error('CHD needs at least 6 columns [iPER L R C hd1 hd2]'); end                


                    if ncol ~= 6 + nAux
                        error('CHD: with %d auxiliary CHD parameters requires 6+%d=%d columns\n',nAux,6+nAux);
                    end

                    %Check if ncol == 7 if any CHDDENSOPT values in column 7 is 1. (swt_v4 manual, p22)
                    if ncol==7
                        if iscell(bcn.BCN)
                            chk7 = any(cellfun(@cols,bcn.BCN)==1);
                        else
                            chk7 = any(bcn.BCN(:,7)==1);
                        end
                        if  any(chk7),
                            fprintf(['You have % columns for CHD, that''s ok if none of the CHDDENSOPT values are 1.\n',...
                                     'In that case you must provide the appropriate density value.\n\n',...
                                     'Lists or sublists with CHDDENSOPT=1, requiring CHDDEN in column: 8'],ncol);

                            error(['CHD needs 6+n colmuns where n is either 0, 1 or 2:\n',...
                                'n=0 if CHDDENSOPT is not set in mf_adapt (i.e. ordinary use of CHD).\n',...
                                'n=1 if CHDDENSOPT is set in mf_adapt like with statement "CHDDENSOPT=1"\n',...
                                '    and the CHDDENSOPT values in column 7 are 0, 2, or 3.\n',...
                                'n=2 if CHDDENSOPT is set and any value in column 7 is 1. Option value 1\n',...
                                '    requires the reference head to be given in column 8.\n',...
                                ' The required BCN list format is (without << and >> of course:)\n',...
                                ' << iPer Layer Row Col  head1 head2 [CHDDENSOPTvalue  [CHDDENrefvalue]] >>']);
                        end
                    end
                end
            end
    end

    bcnp.NPar=0; % parameters not (yet) implemented
    bcn.NPER =NPER;

    MXSS = MXSS + writeBCN(basename,bcn,bcnp);

    clear('bcn'); % and do away with bcn
end

%% ==== THE ASP-file (Doherty's better mf2k MF2KASP) =======================
% TO 100830 to allow using Doherty's MF2KASP, only rewetting options used
asp.SCEN=namSCEN('ASP',nam.PCKG,nam.SCEN);

if asp.SCEN
    % for each simulation
    fprintf('USING MF2KASP, John Doherty''s adapted version of mf2k\n');
    fprintf('Generating asp struct and file\n');
    
    asp.FREE = bas.FREE;

    asp.unit=nam.UNIT(strmatchi('ASP',nam.PCKG));
    asp.ext =nam.EXT {strmatchi('ASP',nam.PCKG)};
    asp.IPESTINT=PESTparvals(strmatchi('IPESTINT',PESTparnams),asp.SCEN);
    asp.INTERP  =PESTparvals(strmatchi('INTERP'  ,PESTparnams),asp.SCEN);
    asp.NOSTOP  =PESTparvals(strmatchi('NOSTOP'  ,PESTparnams),asp.SCEN);
    asp.HDRYBOT =PESTparvals(strmatchi('HDRYBOT' ,PESTparnams),asp.SCEN);
    asp.LIMOP   =PESTparvals(strmatchi('LIMOP'   ,PESTparnams),asp.SCEN);
    asp.MINTHICK=PESTparvals(strmatchi('MINTHICK',PESTparnams),asp.SCEN);
    
    writeASP(basename,asp);
end

%% ==== THE OC-file (output control package) =======================
oc.SCEN=namSCEN('OC',nam.PCKG,nam.SCEN);

if oc.SCEN
    % for each simulation
    fprintf('Generating OC struct and file\n');
    oc.unit   = nam.UNIT(strmatchi('OC',nam.PCKG));
    oc.ext    = nam.EXT {strmatchi('OC',nam.PCKG)};
    oc.IHEDFM = MFLOWparvals(strmatchi('IHEDFM',MFLOWparnams),oc.SCEN);
    oc.IDDNFM = MFLOWparvals(strmatchi('IDDNFM',MFLOWparnams),oc.SCEN);
    oc.IHEDUN = MFLOWparvals(strmatchi('IHEDUN',MFLOWparnams),1); % file unit number for saving heads
    oc.IDDNUN = MFLOWparvals(strmatchi('IDDNUN',MFLOWparnams),1); % file unit number for saving drawdowns
    
    % for each stres period
    % oc.INCODE  =PERvals(:,strmatchi('INCODE',PERnams));  % OBSOLETE as off 130204: always zero
    oc.INCODE  =zeros(size(PERvals(:,1)));                 % now always zero
    oc.IHDDFL = PERvals(:,strmatchi('IHDDFL',PERnams));
    oc.IBUDFL = PERvals(:,strmatchi('IBUDFL',PERnams));
    oc.ICBCFL = PERvals(:,strmatchi('ICBCFL',PERnams));
    
    oc.NPER   = size(PERvals,1);
    oc.NSTP   = PERvals(:,strmatchi('NSTP',PERnams));
    
    oc.Hdpr =PERvals(:,strmatchi('Hdpr',PERnams));
    oc.Ddpr =PERvals(:,strmatchi('Ddpr',PERnams));
    oc.Hdsv =PERvals(:,strmatchi('Hdsv',PERnams));
    oc.Ddsv =PERvals(:,strmatchi('Ddsv',PERnams));

    try
        oc.COMPACT= MFLOWparvals(strmatchi('COMPACT',MFLOWparnams),1);% flag, indicating use of compact budget file or not
        writeOCwords(basename,oc);
    catch ME
        writeOC(basename,oc);
    end

end

%% ==== THE CFP-file (Conduit FLow Process package) =======================
cfp.SCEN=namSCEN('CFP',nam.PCKG,nam.SCEN);

if cfp.SCEN
    cfp.unit=nam.UNIT(strmatchi('CFP',nam.PCKG));
    cfp.ext =nam.EXT {strmatchi('CFP',nam.PCKG)};
    
    cfp.FREE = bas.FREE;

    % Import the file
    [CFPnams ,CFPvals ] = getExcelData(basename,'CFP','V');
    [CFPNnams,CFPNvals] = getExcelData(basename,'CFPNODES','H');
    [CFPPnams,CFPPvals] = getExcelData(basename,'CFPPIPES','H');

    % for each simulation
    cfp.mode       = CFPvals(strmatchi('MODE',CFPnams),1);
    cfp.NNODES     = size(CFPNvals,1);
    cfp.NPIPES     = size(CFPPvals,1);
    cfp.NLAY       = size(LAYparvals,1);
    cfp.temperature= CFPvals(strmatchi('TEMPERATURE',CFPnams),1);
    cfp.LTEMP      = CFPvals(strmatchi('LTEMP',      CFPnams),1);
    cfp.SA_EXCHANGE= CFPvals(strmatchi('SA_EXCHANGE',CFPnams),1);
    cfp.EPSILON    = CFPvals(strmatchi('EPSILON'    ,CFPnams),1);

    % cfp.Connections
    cfp.Nvals=CFPNvals;
    cfp.Pvals=CFPPvals;

    % the neighber node to which this node is connected by a (max to 6) pipe
    NB=zeros(cfp.NNODES,6);
    PB=zeros(cfp.NNODES,6);

    % Connections between pipes and nodes
    pconn=[CFPPvals(:,[1 2 3]); ...
           CFPPvals(:,[1 3 2])...
           ];  % [PipeNr Nd1 Nd2; PipeNr Nd2 Nd1]

    for inode=1:cfp.NNODES
      % Nodes nonnected to this node, we find them by querying the other
      % node of every pipe that is connected to this node
      NNr=CFPNvals(inode,1);

      nb=pconn(pconn(:,2)==NNr,3); %the other end of the pipe connected to this node
      if ~isempty(nb)
        NB(inode,1:numel(nb))=nb;  
      end

      % Pipe connected to this node
      pb = pconn(pconn(:,2)==NNr,1);
      if ~isempty(pb),
          PB(inode,1:numel(pb))=pb;
      end
    end
    cfp.NB=NB;
    cfp.PB=PB;

    % Convergence criteria CFP
    cfp.epsilon =        CFPvals(strmatchi('EPSILON',CFPnams),1);
    cfp.NITER =          CFPvals(strmatchi('NITER',CFPnams),1);
    cfp.RELAX =          CFPvals(strmatchi('RELAX',CFPnams),1);
    cfp.P_NR =           CFPvals(strmatchi('P_NR',CFPnams),1);

    % Pipe properties
    cfp.PIPES =          CFPPvals(:,strmatchi('NO_P',CFPPnams));
    cfp.DIAMETER =       CFPPvals(:,strmatchi('DIAMETER',CFPPnams));
    cfp.TORTUOSITY =     CFPPvals(:,strmatchi('TORTUOSITY',CFPPnams));
    cfp.RHEIGHT =        CFPPvals(:,strmatchi('RHEIGHT',CFPPnams));
    cfp.LCRITREY_P =     CFPPvals(:,strmatchi('LCRITREY_P',CFPPnams));
    cfp.TCRITREY_P =     CFPPvals(:,strmatchi('TCRITREY_P',CFPPnams));

    % Node properties
    cfp.iK_EXCHANGE = strmatchi('K_EXCHANGE',CFPNnams);
    cfp.iN_HEAD     = strmatchi('N_HEAD',CFPNnams);
    cfp.iGEOHEIGHT  = strmatchi('GEOHEIGHT',CFPNnams);

    cfp.NO_N =           CFPNvals(:,strmatchi('NO_N',CFPNnams));
    cfp.K_EXCHANGE =     CFPNvals(:,cfp.iK_EXCHANGE);
    cfp.N_HEAD =         CFPNvals(:,cfp.iN_HEAD);
    cfp.GEOHEIGHT =      CFPNvals(:,cfp.iGEOHEIGHT');

    % Layer properties
    cfp.CL   =           LAYparvals(:,strmatchi('CL',LAYparnams));

    cfp.Lvals =          LAYparvals(:,...
        strmatchi({'VOID','LCRITREY_L','TCRITREY_L'},LAYparnams));

    writeCFP(basename,cfp);

%% ===== THE coc-file (Conduit Flow Process) =====
% Conduit Output Control (COC)
    coc.SCEN=namSCEN('COC',nam.PCKG,nam.SCEN);

    if coc.SCEN && cfp.mode~=2 % if something with pipes
        coc.unit=nam.UNIT(strmatchi('COC',nam.PCKG));
        coc.ext =nam.EXT {strmatchi('COC',nam.PCKG)};

        coc.NNODES       =  cfp.NNODES;
        coc.NODE_NUMBERS =  CFPNvals(:,strmatchi ('NO_N',CFPNnams));
        coc.N_NTS        =  PERvals (strmatchi  ('N_NTS',PERnams));
        coc.NPIPES       =  size(CFPPvals,1);
        coc.PIPE_NUMBERS =  CFPPvals(:,strmatchi('NO_P',CFPPnams));
        coc.T_NTS        =  PERvals(strmatchi ('T_NTS',PERnams));

        writeCOC(basename,coc);
    end
 
% %% ===== THE crch-file (Conduit Flow Process) =====
    crch.SCEN=namSCEN('CRCH',nam.PCKG,nam.SCEN);

    if crch.SCEN && cfp.mode ~=2 % something with pipes
        crch.unit=nam.UNIT(strmatchi('CRCH',nam.PCKG));
        crch.ext =nam.EXT {strmatchi('CRCH',nam.PCKG)};

        crch.IFLAG_CRCH   = CFPvals (strmatchi  ('IFLAG_CRCH',CFPnams));
        crch.NODE_NUMBERS = CFPNvals(:,strmatchi ('NO_N',CFPNnams));
        crch.P_CRCH       = CFPNvals(:,strmatchi ('P_CRCH' , CFPNnams));

        writeCRCH(basename,crch);
    end

end


%% Any solvers on ??
I = ismember(solvers,nam.PCKG);
if sum(I)>1
    error(['%s: You are using more than one solver, namely %s.\n',...
           'REMEDY: Check solvers in NAM sheet'],...
        mfilename,sprintfs(' <<%s>>',solvers(I)));
end
if ~any(I)  && ~any(strmatchi('MT3D',nam.MODEL))
    error(['%s: No solver selected!\n',...
           'This may be due to using a solver that is not supported in the selected model.\n',...
           'For instance, MF2000 and MF2007 does not recognize packages UPW and (solver) NWT.\n',...
           'Therefore, selecting NWT solver with such models does not work.\n'....
           'REMEDY: in the NAM sheet select a solver that is recognized for the\n',...
           '        selected model. PCG, SIP, SOR and DE4 generally work.\n',...
           '        PCGN may be an alternative for non linear problems in MF2K and MF2005\n',...
           '        NWT only works with mf2005nwt (together with UPW).'],mfilename);
end

%% ===== THE pcg-file (P Conjugate solver Package) =====
pcg.SCEN=namSCEN('PCG',nam.PCKG,nam.SCEN);

if pcg.SCEN    
    fprintf('Generating PCG struct and file\n');
        
    pcg.unit=nam.UNIT(strmatchi('PCG',nam.PCKG,'exact'));
    pcg.ext =nam.EXT {strmatchi('PCG',nam.PCKG,'exact')};
    
    fetch = @(nm) MFLOWparvals(strmatchi(nm,MFLOWparnams,'exact'),1);
    
    variables = {'MXITER','ITER1','NPCOND','HCLOSE','RCLOSE', ...
                 'RELAX','NBPOL','IPRPCG','MUTPCG','DAMP'};

    for i=1:numel(variables)
        pcg.(variables{i}) = fetch(variables{i});
    end

    pcg.RELAX = pcg.RELAX(1); % conflicts with pcgn
    pcg.DAMP  = pcg.DAMP(1);  % conflicts with pcgn
    
    writePCG(basename,pcg);
end

%% ===== THE pcgn-file (P Conjugate solver Package with nonlinear improvements) =====
pcgn.SCEN=namSCEN('PCGN',nam.PCKG,nam.SCEN);

if pcgn.SCEN
    fprintf('Generating PCG struct and file\n');
    
    pcgn.unit=nam.UNIT(strmatchi('PCGN',nam.PCKG,'exact'));
    pcgn.ext =nam.EXT {strmatchi('PCGN',nam.PCKG,'exact')};

    fetch = @(nm) MFLOWparvals(strmatchi(nm,MFLOWparnams,'exact'),1);

    variables = {'ITER_MO','ITER_MI','CLOSE_R','CLOSE_H',...
                  'RELAX','IFILL','UNIT_PC','UNIT_TS',...
                  'ADAMP','DAMP','DAMP_LB','RARE_D','CHGLIMIT',...
                  'ACNVG','CNVG_LB','MCNVG','RATE_C','IPUNIT'};
    
    for i = 1:numel(variables)
        pcgn.(variables{i}) = fetch(variables{i});
    end

    pcg.RELAX = pcg.RELAX(end); % conflicts with pcg
    pcg.DAMP  = pcg.DAMP(end);  % conflicts with pcg

    writePCGN(basename,pcgn);
end


%% ===== THE nwt-file (NWT package (Newton Solver Package)
nwt.SCEN=namSCEN('NWT',nam.PCKG,nam.SCEN);

if nwt.SCEN
    fprintf('Generating NWT struct and file\n');
    if ~namSCEN('UPW',nam.PCKG,nam.SCEN)
        error(['%s: You must use NWT together with UPW.\n',...
               'REMEDY: Check NAM sheet.\n',...
               '        Use other solver if you don''t use UPW\n',...
               '        USE NWT sovler if you use UPW.'],mfilname);
    end
    nwt.unit=nam.UNIT(strmatchi('NWT',nam.PCKG,'exact'));
    nwt.ext =nam.EXT {strmatchi('NWT',nam.PCKG,'exact')};
    
    nwt.HEADTOL    = MFLOWparvals(strmatchi('HEADTOL'   ,MFLOWparnams,'exact'),1);
    nwt.FLUXTOL    = MFLOWparvals(strmatchi('FLUXTOL'   ,MFLOWparnams,'exact'),1);
    nwt.MAXITEROUT = MFLOWparvals(strmatchi('MAXITEROUT',MFLOWparnams,'exact'),1);
    nwt.THICKFACT  = MFLOWparvals(strmatchi('THICKFACT' ,MFLOWparnams,'exact'),1);
    nwt.LINMETH    = MFLOWparvals(strmatchi('LINMETH'   ,MFLOWparnams,'exact'),1);
    nwt.IPRNWT     = MFLOWparvals(strmatchi('IPRNWT'    ,MFLOWparnams,'exact'),1);
    nwt.IBOTAV     = MFLOWparvals(strmatchi('IBOTAV'    ,MFLOWparnams,'exact'),1);
    nwt.OPTIONS    = MFLOWparvals(strmatchi('OPTIONS'   ,MFLOWparnams,'exact'),1);

    if nwt.OPTIONS == 4
    
        nwt.DBDTHETA    = MFLOWparvals(strmatchi('DBDTHETA'   ,MFLOWparnams,'exact'),1);
        nwt.DBDKAPPA    = MFLOWparvals(strmatchi('DBDKAPPA'   ,MFLOWparnams,'exact'),1);
        nwt.DBDGAMMA    = MFLOWparvals(strmatchi('DBDGAMMA'   ,MFLOWparnams,'exact'),1);
        nwt.MOMFACT     = MFLOWparvals(strmatchi('MOMFACT'    ,MFLOWparnams,'exact'),1);
        nwt.BACKFLAG    = MFLOWparvals(strmatchi('BACKFLAG'   ,MFLOWparnams,'exact'),1);
        nwt.MAXBACKITER = MFLOWparvals(strmatchi('MAXBACKITER',MFLOWparnams,'exact'),1);
        nwt.BACKTOL     = MFLOWparvals(strmatchi('BACKTOL'    ,MFLOWparnams,'exact'),1);
        nwt.BACKREDUCE  = MFLOWparvals(strmatchi('BACKREDUCE' ,MFLOWparnams,'exact'),1);

        switch nwt.LINMETH
           case 1
               nwt.MAXITINNER = MFLOWparvals(strmatchi('MAXITINNER',MFLOWparnams,'exact'),1);
               nwt.ILUMETHOD  = MFLOWparvals(strmatchi('ILUMETHOD' ,MFLOWparnams,'exact'),1);
               nwt.LEVFILL    = MFLOWparvals(strmatchi('LEVFILL'   ,MFLOWparnams,'exact'),1);
               nwt.STOPTOL    = MFLOWparvals(strmatchi('STOPTOL'   ,MFLOWparnams,'exact'),1);
               nwt.MSDR       = MFLOWparvals(strmatchi('MSDR'      ,MFLOWparnams,'exact'),1);
           case 2
               nwt.IACL       = MFLOWparvals(strmatchi('IACL'      ,MFLOWparnams,'exact'),1);
               nwt.NORDER     = MFLOWparvals(strmatchi('NORDER'    ,MFLOWparnams,'exact'),1);
               nwt.LEVEL      = MFLOWparvals(strmatchi('LEVEL'     ,MFLOWparnams,'exact'),1);
               nwt.NORTH      = MFLOWparvals(strmatchi('NORTH'     ,MFLOWparnams,'exact'),1);
               nwt.IREDSYS    = MFLOWparvals(strmatchi('IREDSYS'   ,MFLOWparnams,'exact'),1);
               nwt.RRCTOLS    = MFLOWparvals(strmatchi('RRCTOLS'   ,MFLOWparnams,'exact'),1);
               nwt.IDROPTOL   = MFLOWparvals(strmatchi('IDROPTOL'  ,MFLOWparnams,'exact'),1);
               nwt.EPSRN      = MFLOWparvals(strmatchi('EPSRN'     ,MFLOWparnams,'exact'),1);
               nwt.HCLOSEXMD  = MFLOWparvals(strmatchi('HCLOSEXMD' ,MFLOWparnams,'exact'),1);
               nwt.MXITERXMD  = MFLOWparvals(strmatchi('MXITERXMD' ,MFLOWparnams,'exact'),1);
           otherwise
               error('%s: illegal LINMETH = %d, must be 1 or 2',mfilename,nwt.LINMETH);
        end
    end
    writeNWT(basename,nwt);
end


%% ===== THE de4-file (P Conjugate solver Package) =====
de4.SCEN=namSCEN('DE4',nam.PCKG,nam.SCEN);

if de4.SCEN
    fprintf('Generating PCG struct and file\n');
    de4.unit  =nam.UNIT(    strmatchi('DE4'   ,nam.PCKG    ,'exact'));
    de4.ext   =nam.EXT {    strmatchi('DE4'   ,nam.PCKG    ,'exact')};
    de4.ITMX  =MFLOWparvals(strmatchi('ITMX'  ,MFLOWparnams,'exact'),de4.SCEN);
    de4.MXUP  =MFLOWparvals(strmatchi('MXUP'  ,MFLOWparnams,'exact'),de4.SCEN);
    de4.MXLOW =MFLOWparvals(strmatchi('MXLOW' ,MFLOWparnams,'exact'),de4.SCEN);
    de4.MXBW  =MFLOWparvals(strmatchi('MXBW'  ,MFLOWparnams,'exact'),de4.SCEN);
    de4.IFREQ =MFLOWparvals(strmatchi('IFREQ' ,MFLOWparnams,'exact'),de4.SCEN);
    de4.MUTD4 =MFLOWparvals(strmatchi('MUTD4' ,MFLOWparnams,'exact'),de4.SCEN);
    de4.ACCL  =MFLOWparvals(strmatchi('ACCL'  ,MFLOWparnams,'exact'),de4.SCEN);
    de4.HCLOSE=MFLOWparvals(strmatchi('HCLOSE',MFLOWparnams,'exact'),de4.SCEN);
    de4.IPRD4 =MFLOWparvals(strmatchi('IPRD4' ,MFLOWparnams,'exact'),de4.SCEN);
    
    writeDE4(basename,de4);
end

%% ===== THE SIP-file (Strongly Implicit Procedure Package) =====
sip.SCEN=namSCEN('SIP',nam.PCKG,nam.SCEN);

if sip.SCEN
    fprintf('Generating SIP struct and file\n');
    sip.unit   =nam.UNIT(    strmatchi('SIP'    ,nam.PCKG    ,'exact'));
    sip.ext    =nam.EXT {    strmatchi('SIP'    ,nam.PCKG    ,'exact')};
    sip.MXITER =MFLOWparvals(strmatchi('MXITER' ,MFLOWparnams,'exact'),sip.SCEN);
    sip.NPARM  =MFLOWparvals(strmatchi('NPARM'  ,MFLOWparnams,'exact'),sip.SCEN);
    sip.ACCL   =MFLOWparvals(strmatchi('ACCL'   ,MFLOWparnams,'exact'),sip.SCEN);
    sip.HCLOSE =MFLOWparvals(strmatchi('HCLOSE' ,MFLOWparnams,'exact'),sip.SCEN);
    sip.IPCALC =MFLOWparvals(strmatchi('IPCALC' ,MFLOWparnams,'exact'),sip.SCEN);
    sip.WSEED  =MFLOWparvals(strmatchi('WSEED'  ,MFLOWparnams,'exact'),sip.SCEN);
    sip.IPRSIP =MFLOWparvals(strmatchi('IPRSIP' ,MFLOWparnams,'exact'),sip.SCEN);
    
    writeSIP(basename,sip);
end

%% ===== THE SOR-file (Sliced Successive Over Relaxation Package) =====
sor.SCEN=namSCEN('SOR',nam.PCKG,nam.SCEN);

if sor.SCEN
    fprintf('Generating SOR struct and file\n');
    sor.unit   =nam.UNIT(    strmatchi('SOR'    ,nam.PCKG    ,'exact'));
    sor.ext    =nam.EXT {    strmatchi('SOR'    ,nam.PCKG    ,'exact')};
    sor.MXITER =MFLOWparvals(strmatchi('MXITER' ,MFLOWparnams,'exact'),sor.SCEN);
    sor.ACCL   =MFLOWparvals(strmatchi('ACCL'   ,MFLOWparnams,'exact'),sor.SCEN);
    sor.HCLOSE =MFLOWparvals(strmatchi('HCLOSE' ,MFLOWparnams,'exact'),sor.SCEN);
    sor.IPRSOR =MFLOWparvals(strmatchi('IPRSOR' ,MFLOWparnams,'exact'),sor.SCEN);
    
    writeSOR(basename,sor);
end

%% ===== THE LMT6 file (Link file generation LMT6) =====
lmt.SCEN=namSCEN('LMT',nam.PCKG,nam.SCEN);

if lmt.SCEN
    lmt.unit=nam.UNIT(strmatchi('LMT6',nam.PCKG,'exact'));  % must not be used !!
    lmt.ext =nam.EXT {strmatchi('LMT6',nam.PCKG,'exact')};
    lmt.OUTPUT_FILE_NAME   = [basename '.FTL'];
    lmt.OUTPUT_FILE_HEADER ='EXTENDED';
    lmt.OUTPUT_FILE_FORMAT ='UNFORMATTED';

    writeLMT(basename,lmt);
end

%% ===== END OF INPUT FOR THE FLOW PROCESs ==================

%% ===== MODPATH INPUT=================================
% mfLab implements MODPATH version 6.0
% David W. Pollock (2012)
% User Guide for MODPATH Version 6 ? A Particle-Tracking Model for MODFLOW
% Chapter 41 of Section A, Groundwater Book 6, Modeling Techniques
% Older versions are not supported.

%% ===== THE modpath v 6.0 simlation file

%% End mdpath v 6.0 simulation file
I = strmatchi('MODPATH',nam.MODEL); I=I(1);

if  I(1)
    fprintf('Generating MODPATH v6.0 simulation struct and file\n');

    [~,mpath_basename] = fileparts(nam.mdlpath{I(1)});

    %% Data from worksheet MPATH
    [mpathNams,mpathVals,mpathTxtHdr,mpathTxtVals] = getExcelData(basename,'MPATH','v');

    %% File names
    pth.dims             = GRID.size;
    pth.LAYCBD           = GRID.LAYCBD;
    pth.simulationFile   = mpathTxtVals{strmatchi('simulationFile',mpathTxtHdr),1};
    pth.nameFile         = mpathTxtVals{strmatchi('nameFile',mpathTxtHdr),1};
    pth.listingFile      = mpathTxtVals{strmatchi('listingFile',mpathTxtHdr),1};
    pth.startingLocationsFile = mpathTxtVals{strmatchi('startingLocationsFile',mpathTxtHdr),1};
    pth.endPointsFile    = mpathTxtVals{strmatchi('endPointsFile',mpathTxtHdr),1};
    pth.pathLinesFile    = mpathTxtVals{strmatchi('pathLinesFile',mpathTxtHdr),1};
    pth.timeSeriesFile   = mpathTxtVals{strmatchi('timeSeriesFile',mpathTxtHdr),1};
    pth.observationsFile = mpathTxtVals{strmatchi('observationsFile',mpathTxtHdr),1};
    pth.traceFile        = mpathTxtVals{strmatchi('traceFile',mpathTxtHdr),1};

    %% Options
    %3
    pth.simulationType   = mpathVals(strmatchi('simulationType',mpathNams),1);
    
    %% TrackingDirection
    pth.trackingDirection  = mpathVals(strmatchi('TrackingDirection',mpathNams),1);
    switch pth.trackingDirection
        case 1 % foward
        case 2 % backward
        otherwise
            error('trackingDirection must be 1(forward) or 2(backward), see worksheet MPATH');
    end
    
    %% weakSinkOption
    pth.weakSinkOption   =  mpathVals(strmatchi('WeakSinkOption',mpathNams,1));
    switch pth.weakSinkOption
        case 1
        case 2
        otherwise
            error('Unknown option <<%d>> for weakSinkOption, use 1..2, see worksheet MPATH',...
                mpth.weakSinkOption);
    end
    
    %% weakSourceOption
    pth.weakSourceOption =  mpathVals(strmatchi('WeakSourceOption',mpathNams),1);
    switch pth.weakSinkOption
        case 1
        case 2
        otherwise
            error('Unknown option <<%d>> for weakSourceOption, use 1..2, see worksheet MPATH',...
                mpth.weakSinkOption);
    end
    
    %% ReferenceTimeOption
    pth.referenceTimeOption  =  mpathVals(strmatchi('referenceTimeOption',mpathNams),1);
    I = strmatchi('referenceTime',mpathNams,'exact');
    switch pth.referenceTimeOption
        case 1
            pth.referenceTime = mpathVals(I(1),1);
            if pth.referenceTime<0
                error('referenceTime is %g, must be >=0, see worksheet MPATH',...
                    pth.referenceTime);
            end
        case 2
            pth.period       = mpathVals(I(2),1);
            if pth.period<1
                error('period nr<1, see worksheet MPATH');
            end
            pth.step         = mpathVals(I(2),2);
            if pth.step<1
                error('step nr<1, see worksheet MPATH');
            end
            pth.timeFraction = mpathVals(I(2),3);
            if pth.timeFraction<0 || pth.timeFraction>1
                eror('timeFraction must be between 0 and 1, see worksheet MPATH');
            end
            pth.referenceTime= -1;
        otherwise
    end
    
    %% StopOption 
    pth.stopOption = mpathVals(strmatchi('stopOption',mpathNams),1);
    switch pth.stopOption
        case 1
            % Stop at end or beginning of Modflow simulation.
        case 2
            % Extend initial or final steady state Modflow time step
            % as far as necessary to track all particles to their termination points.
        case 3
            % Specify stop time
            pth.stopTime = mpathVals(strmatchi('stopTime',mpathNams),1);
        otherwise
            error('Unknown option <<%d>> for stopOption (see sheet MPATH)',pth.stopOption);
    end
    
    %% ParticleGenerationOption
    pth.particleGenerationOption = mpathVals(strmatchi('particleGenerationOption',mpathNams),1);
    switch pth.particleGenerationOption
        case 1  % Put particle generation info in simulation file
        case 2  % Put particle generation info in particleStartingLocations file
        otherwise
            error('Unknown value <<%d>> for particleGenerationOption, use 1..2, see worksheet MPATH');
    end
    
    %% TimePointOption
    pth.timePointOption = mpathVals(strmatchi('timePointOption',mpathNams),1);
    pth.timePointCount  = mpathVals(strmatchi('timePointCount',mpathNams),1);
    if pth.timePointCount<1
        error('timePointCount<1, see worksheet MPATH');
    end
    switch pth.timePointOption
        case 1
        case 2
            pth.releaseTimeIncrement = mpathVals(strmatchi('releaseTimeIncrement',mpathNams),1);
            if pth.releaseTimeIncrement<=0
                error('releaseTimeIncremetn<=0, see worksheet MPATH');
            end
        case 3
            pth.timePoints = mpathVals(strmatchi('timePoints',mpathNams),:);
            if pth.timePointCount> numel(pth.timePoints)
                error('timePointCount<number of available timePoints. See worksheet MPATH');
            end
            pth.timePoints = pth.timePoints(1:pth.timePointCount);
        otherwise
            error('Illegal timePointOption <<%d>> must be 1..3, see worksheet MPATH',...
                pth.timePointOption);
    end
        
    %% BudgetOutputOption
    pth.budgetOutputOption = mpathVals(strmatchi('BudgetOutputOption',mpathNams),1);
    switch pth.budgetOutputOption
        case 1
        case 2
            % print budget summar of all cells in Listing file
        case 3
            try
                pth.budgetOutputCells = BUDGETCELLS;
            catch ME
                fprintf(1,ME.message); fprintf('\n');
                error(['For budgetOutputOption ==3, you must specify a list of\n',...
                    'output cells in the workspace with the name BUDGETOUPUTCELLS']);
            end
        case 4
            % Trace is active
        otherwise
            error('unknown option <<%s>> for budgetOutputOption, use 1..4, see sheet MPATH',...
                pth.budgetOutputOption);
    end
    
    %% zoneArrayOption (to define a stopzone)
    pth.zoneArrayOption = mpathVals(strmatchi('zoneArrayOption',mpathNams),1);
    switch pth.zoneArrayOption
        case 1
            % no zone array will be read
        case 2
            try
                pth.ZONE = ZONE;
            catch ME
                fprintf(1,ME.message); fprintf('\n');
                error(['with zoneArrayOption == %d, you must specify a 3D arrray ZONE\n',...
                    'in the workshpace to define the stop zone for Modpath'],...
                    pth.zoneArrayOption);
            end
            % from spreadsheet
            pth.stopZone = mpathVals(strmatchi('stopZone',mpathNams),1);
        otherwise
            error('unknown value <<%d>> for zoneArrayOption, use 1..2, see worksheet MPATH',...
                pth.zoneArrayOption);            
    end
    
    %% retardationOption
    
    pth.retardationOption = mpathVals(strmatchi('retardationOption',mpathNams),1);
    switch pth.retardationOption
        case 1 % no retardation
        case 2
            try
                pth.RETARDATION = RETARDATION;
            catch ME
                fprintf(1,ME.message); fprintf('\n');
                error(['If retardationOption==2, you must specify\n',...
                    'a 3D array RETARDATION in the workspace, see worksheet MPATH']);
            end
            if any(GRID.LAYCBD)
                try
                    pth.RETARDATIONCB = RETARDATIONCB;
                catch ME
                    fprintf(1,ME.message); fprintf('\n');
                    error(['if retardationOption==2 and any(gridObj.LAYCBD>0), you must specify\n'...
                        ' a 3D array RETARDATIONCB in the workspace, see worksheet MPATH']);
                end
            end
        otherwise
    end

    %% AdvectiveObservationsOption
    pth.advectiveObservationsOption = ...
        mpathVals(strmatchi('advectiveObservationsOption',mpathNams),1);
    switch pth.advectiveObservationsOption
        case 1  % no advective obserations are computed or saved
        case 2  % adv. obs. saved for all time points
        case 3  % adv. obs. saved for final time point only.
        otherwise
            error(['Unknown option <<%d>> for advectiveObservationsOption,\n',...
                'Use 1..3, see workshet MPATH'],pth.advectiveObservationsOption);
    end
    
    %% add mpath_particleGroupObj's to pth
    eval(['pth.particleGroup = ' particleGroupNm]);
    
    %% Further implied by input of mfLab and not specified here
    % [minLayer minRow minColumn, maxLayer maxRow maxColumn]
    % mask(NCOL,NROW)
    % maskLayer
    % Facecount
    % IFace
    % ParticleLayerCount ParticleRowCount ParticleColumnCount
    
    %%
    writeMPTH6(mpath_basename,pth);
       
    %% MPBAS MODPATH basic Data file data
    
    I = strmatchi('MPBAS',nam.PCKG);
    mpbas.unit = nam.UNIT(I(1));
    mpbas.ext  = nam.EXT{I(1)};
    
    mpbas.HNOFLO = bas.HNOFLO;
    try
        mpbas.HDRY   = bas.HDRY;
    catch ME
        mpbas.HDRY   = bas.HNOFLO;
    end
    
    mpbas.NLAY   = GRID.Nlay;

    % Item 2 DefaultIFaceCount
    mpbas.budgetLabel  = mpathTxtVals(strmatchi('budgetLabel' ,mpathTxtHdr),:);
    mpbas.defaultIFace = mpathVals(   strmatchi('defaultIFace',mpathNams),:);
    
    % item5
    try
        mpbas.LAYTYP=lpf.LAYCON;
    catch ME
        mpbas.LAYTYP=bcf.LAYCON;
    end
    
    % item6
    mpbas.IBOUND = bas.IBOUND;

    %9 Effective porosity for the layers
    mpbas.porosity   = PEFF;

    % if any LAYCBD, then PORCB must be specified for Modpath
    mpbas.LAYCBD = GRID.LAYCBD;
    
    if any(GRID.LAYCBD)
        try
            mpbas.porosityCB = PORCB;
        catch ME
            fprintf(2,'Can''t fined PORCB, effective porosity for confining beds\n');
            throw(ME);
        end
    end
    writeMPBAS(basename,mpbas);

end % of modpath

%% ==== MT3DMS input ==================================
% partly from the MT3D sheet in the spreadsheet


%% THE BTN-file (Basic Transprot Process for MT3D)
btn.SCEN=namSCEN('BTN',nam.PCKG,nam.SCEN);

if btn.SCEN
    btn.FREE  = FREE;
    btn.NCOMP = MT3Dvals(strmatchi('NCOMP',MT3Dnams),btn.SCEN);
    btn.MCOMP = MT3Dvals(strmatchi('MCOMP',MT3Dnams),btn.SCEN);
    
    fprintf('Generating Basic Transport Process struct\n');
    btn.unit = nam.UNIT(strmatchi('BTN',nam.PCKG));
    btn.ext  = nam.EXT {strmatchi('BTN',nam.PCKG)};
    %3
    btn.NLAY = GRID.Nlay;   % # layers
    btn.NROW = GRID.Ny;   % # rows
    btn.NCOL = GRID.Nx;   % # columns
    btn.NPER = NPER;   % # stress periods
    
    %4
    TUNIT = MFLOWparvals(strmatchi('ITMUNI',MFLOWparnams),1);
    LUNIT = MFLOWparvals(strmatchi('LENUNI',MFLOWparnams),1);
    switch TUNIT
        case 1, btn.TUNIT = 'SEC';
        case 2, btn.TUNUT = 'MIN';
        case 3, btn.TUNIT = 'HOUR';
        case 4, btn.TUNIT = 'DAYS';
        case 5, btn.TUNIT = 'YEAR';
        otherwise
            btn.TUNIT='????';
    end
    switch LUNIT
        case 1, btn.LUNIT = 'FEET';
        case 2, btn.LUNIT = 'M';
        case 3, btn.LUNIT = 'CM';
        otherwise
            btn.LUNIT = '????';
    end
    btn.MUNIT='KG';
    
    %5
    %Packge-use flags: (ADV DSP SSM RCT GCG XXX XXX XXX XXX XXX)';
    adv.SCEN=strmatchi('ADV',nam.PCKG);
    dsp.SCEN=strmatchi('DSP',nam.PCKG);
    ssm.SCEN=strmatchi('SSM',nam.PCKG);
    rct.SCEN=strmatchi('RCT',nam.PCKG);
    gcg.SCEN=strmatchi('GCG',nam.PCKG);

    % active flags then become
    btn.TRNOP=[adv.SCEN,dsp.SCEN,ssm.SCEN,rct.SCEN,gcg.SCEN,0,0,0,0,0];
    
    %6
    btn.LAYCON=LAYparvals(:,strmatchi('LAYCON',LAYparnams));
    %7
    btn.DELR=GRID.dx;
    %8
    btn.DELC=GRID.dy;
    %9 %10
    btn.Z=GRID.Z;
    %11
    if GRID.AXIAL, btn.PRSITY=GRID.TWOPIR.*PEFF; else btn.PRSITY=PEFF; end
    %12
    btn.ICBUND = ICBUND;  % 0 inactive, 1 active, -1 fixed concentration 
    %13
    btn.STCONC=STCONC;
    %14
    btn.CINACT = MT3Dvals(strmatchi('CINACT',MT3Dnams),btn.SCEN);
    btn.THKMIN = MT3Dvals(strmatchi('THKMIN',MT3Dnams),btn.SCEN);
    %15
    btn.IFMTCN = MT3Dvals(strmatchi('IFMTCN',MT3Dnams),btn.SCEN);
    btn.IFMTNP = MT3Dvals(strmatchi('IFMTNP',MT3Dnams),btn.SCEN);
    btn.IFMTRF = MT3Dvals(strmatchi('IFMTRF',MT3Dnams),btn.SCEN);
    btn.IFMTDP = MT3Dvals(strmatchi('IFMTDP',MT3Dnams),btn.SCEN);
    btn.SAVUCN = MT3Dvals(strmatchi('SAVUCN',MT3Dnams),btn.SCEN);
    %16 %17
    % Note if TIMPRS(1)==0 it will be truncated in writeBTN
    btn.NPRS   = MT3Dvals(strmatchi('NPRS',MT3Dnams),btn.SCEN);
    btn.TIMPRS = MT3Dvals(strmatchi('TIMPRS',MT3Dnams),:); % record with times to print output
    
    % The MT3D option btn.NPRS is quite useless normally if one wants to
    % synchronize output of MT3DMS with that of the heads and the
    % cell-by-cell flow terms written to the budget file. In mfLab, if
    % btn.NPRS>0, we synchronize with the head output.
    if btn.NPRS>0
        btn.TIMPRS=MT3D_timprs(PERnams,PERvals);
        btn.NPRS  =numel(btn.TIMPRS);
    end
    
    %18
    btn.NPROBS=MT3Dvals(strmatchi('NPROBS',MT3Dnams),btn.SCEN);
     %19  skipped, are the locations of the observation points
    %% OBSERVATION POINTS FOR CONCENTRATION
    %if btn.NOBS>0
    try
        [~,btn.OBS]=getExcelData(XLSF,'BTNOBS','Horizontal');
        if size(btn.OBS,1)>0
            fprintf('%d conc observation points read from worksheet BTNOBS\n',size(btn.OBS,1));
        end
    catch ME
        btn.OBS=[];
        fprintf('No worksheet with concentration observation points.\n');
    end
    if exist('BTNOBS','var')
        fprintf('%d conc observation points obtained from mf_adapt\n',size(BTNOBS,1));
        btn.OBS=[btn.OBS;BTNOBS];
    end
    %end
    %20
    btn.CHKMAS = MT3Dvals(strmatchi('CHKMAS',MT3Dnams),btn.SCEN);
    btn.NPRMAS = MT3Dvals(strmatchi('NPRMAS',MT3Dnams),btn.SCEN);
    %21
    btn.PERLEN = PERvals(:,strmatchi('PERLEN',PERnams));
    btn.NSTP   = PERvals(:,strmatchi('NSTP',  PERnams));
    btn.TSMULT = PERvals(:,strmatchi('TSMULT',PERnams));
    %23
    btn.DT0    = PERvals(:,strmatchi('DT0'    ,PERnams)); % i.e. let the model choose
    btn.MXSTRN = PERvals(:,strmatchi('MXSTRN' ,PERnams)); % maximum number of transport steps within one flow time step
    btn.TTSMULT= PERvals(:,strmatchi('TTSMULT',PERnams)); % succuessive transport time step multiplier
    btn.TTSMAX = PERvals(:,strmatchi('TTSMAX' ,PERnams)); %  use default, set no maximum
    
    writeBTN(basename,btn);
end


%% ===== THE ADV-file (Advection process) ==========
adv.SCEN=namSCEN('ADV',nam.PCKG,nam.SCEN);

if adv.SCEN   
    fprintf('Generating Advection Process struct\n');
    adv.FREE=FREE;
    adv.unit=nam.UNIT(strmatchi('ADV',nam.PCKG));
    adv.ext =nam.EXT {strmatchi('ADV',nam.PCKG)};
    %3
    adv.MIXELM=MT3Dvals(strmatchi('MIXELM' ,MT3Dnams),adv.SCEN);   % standard finite difference method
    adv.PERCEL=MT3Dvals(strmatchi('PERCEL' ,MT3Dnams),adv.SCEN);
    adv.MXPART=MT3Dvals(strmatchi('MXPART' ,MT3Dnams),adv.SCEN);
    adv.NADVFD=MT3Dvals(strmatchi('NADVFD' ,MT3Dnams),adv.SCEN);
    %2
    adv.ITRACK=MT3Dvals(strmatchi('ITRACK' ,MT3Dnams),adv.SCEN);   % particle tracking algorith flag 3=good compromise
    adv.WD    =MT3Dvals(strmatchi('WD'     ,MT3Dnams),adv.SCEN);   % concentration weighting factor, 05 should be adequat
    %3
    adv.DCEPS =MT3Dvals(strmatchi('DCEPS'  ,MT3Dnams),adv.SCEN);
    adv.NPLANE=MT3Dvals(strmatchi('NPLANE' ,MT3Dnams),adv.SCEN);
    adv.NPL   =MT3Dvals(strmatchi('NPL'    ,MT3Dnams,'exact'),adv.SCEN);
    adv.NPH   =MT3Dvals(strmatchi('NPH'    ,MT3Dnams),adv.SCEN);
    adv.NPMIN =MT3Dvals(strmatchi('NPMIN'  ,MT3Dnams),adv.SCEN);
    adv.NPMAX =MT3Dvals(strmatchi('NPMAX'  ,MT3Dnams),adv.SCEN);
    %4
    adv.INTERP =MT3Dvals(strmatchi('INTERP',MT3Dnams),adv.SCEN);
    adv.NLSINK =MT3Dvals(strmatchi('NLSINK',MT3Dnams),adv.SCEN);
    adv.NPSINK =MT3Dvals(strmatchi('NPSINK',MT3Dnams),adv.SCEN);
    
    adv.DCHMOC =MT3Dvals(strmatchi('DCHMOC',MT3Dnams),adv.SCEN);
    writeADV(basename,adv);
end


%% ===== THE DSP-file (Dispersion process package) =============
dsp.SCEN=namSCEN('DSP',nam.PCKG,nam.SCEN);

if dsp.SCEN
   fprintf('Generating Dispersion Process struct\n');
    dsp.unit=nam.UNIT(strmatchi('DSP',nam.PCKG));
    dsp.ext =nam.EXT {strmatchi('DSP',nam.PCKG)};
    dsp.NCOMP=btn.NCOMP;
    %3
    dsp.NLAY=GRID.Nlay;
    dsp.NROW=GRID.Ny;
    dsp.NCOL=GRID.Nx;
    
    dsp.MultiDiffusion=MT3Dvals(strmatchi('MULTIDIFFUSION',MT3Dnams),dsp.SCEN);
    
    % if AL exists in mf_adapt, use it. It must be a cell array iLay long
    % each cell has the data for one layer (NROW,NCOL) or a single value
    if exist('AL','var'),
        dsp.AL=AL;
    else
        dsp.AL=cell(dsp.NLAY,1);
    
        dsp.al    =LAYparvals(:,strmatchi('AL'    ,LAYparnams));
        for iLay=1:dsp.NLAY
            if dsp.al(iLay)>=0,  % use spreadsheet value if value >=0
                dsp.AL{iLay}=dsp.al(iLay);
            end
        end
    end

    dsp.TRPT  =LAYparvals(:,strmatchi('TRPT'  ,LAYparnams));  % aTH/aL
    dsp.TRPV  =LAYparvals(:,strmatchi('TRPV'  ,LAYparnams));  % aTV/aL
        dmcoef=LAYparvals(:,strmatchi('DMCOEF',LAYparnams));  % effctief molecular diffusion coefficient
    if size(dmcoef,2)>dsp.NCOMP
        dmcoef=dmcoef(:,1:dsp.NCOMP);
    end
    % if DMCOEF exists in the workspace, we will use it instead of the
    % values in the LAY worksheet.
    if exist('DMCOEF','var'),
        %In principle this assumes DMCOEF is a cell array with one full 3D
        %array per species (==cell)
        %not being stringent about how DMDOEF is specified provides some
        %flexibility because warray in writeDSP will work for a single
        %value per layer as for a full 3D array.
        if iscell(DMCOEF)
            dsp.DMCOEF=DMCOEF;
        else
            % if not a cell array then make it a cell array
            dsp.DMCOEF = {DMCOEF};
        end
    else
        % if not, we get the values from the layworkhsheet. This always
        % yields one column per species while the number of layers is fixed.
        % We will put each column pertaining to a different species in its
        % own cell. Each cell therefore, contains a column of DMCOEF
        % values, that correspond to the layers.
        dsp.DMCOEF=cell(1,dsp.NCOMP);
        for iComp=1:dsp.NCOMP
            dsp.DMCOEF{1,iComp}=XS(dmcoef(:,iComp));
        end
    end
    
    writeDSP(basename,dsp);
end


%% ===== THE SSM-file (Source-sink mixing process package) =====
ssm.SCEN=namSCEN('SSM',nam.PCKG,nam.SCEN);

if ssm.SCEN
    fprintf('Generating Source-Sink Mixing Process struct\n');
    ssm.FREE     = FREE;
    ssm.basename = basename;
    ssm.unit=nam.UNIT(strmatchi('SSM',nam.PCKG));
    ssm.ext =nam.EXT {strmatchi('SSM',nam.PCKG)};
    
    ssm.NPER=NPER;
    ssm.NROW=GRID.Ny;
    ssm.NCOL=GRID.Nx;
    ssm.NLAY=GRID.Nlay;
    ssm.NCOMP=btn.NCOMP;
    ssm.MCOMP=btn.MCOMP;  % don't know if needed here ## check 
    ssm.FRCH = strmatchi('RCH',nam.PCKG,'exact');
    ssm.FEVT = strmatchi('EVT',nam.PCKG,'exact') || strmatchi('ETS',nam.PCKG,'exact');

    %D1 FWEL FDRN FRCH FEVT FRIV FGHB FCHD FMNW (FNEW(n), n=1:3) 10I2
    if ssm.FRCH % recharge concentration has to be defined
        ssm.INCRCH =PERvals(:,strmatchi('INCRCH',PERnams));
        if any(ssm.INCRCH==0), % Then we need to define CRCH in mf_adapt
            try
                ssm.CRCH=CRCH;  % has CRCH been defined in mf_adapt?
            catch ME
               error('If any IN%s in the PER worksheet = 0 then you must define CRCH{NLAY,NCOMP} in mf_adapt','CRCH');
            end
            for iComp=1:ssm.NCOMP
                crch=PERvals(:,strmatchi(sprintf('CRCH_%d',iComp),PERnams));
                for iPer=1:ssm.NPER
                  if ssm.INCRCH(iPer)>0, ssm.CRCH{iPer,iComp}=crch(iPer); end
                end
            end
        else
            ssm.CRCH=cell(ssm.NPER,ssm.NCOMP);
            for iComp=1:ssm.NCOMP
                crch=PERvals(:,strmatchi(sprintf('CRCH_%d',iComp),PERnams));
                for iPer=1:ssm.NPER
                  ssm.CRCH{iPer,iComp}=crch(iPer);
                end
            end
        end
    end
    
    if ssm.FEVT
        ssm.INCEVT =PERvals(:,strmatchi('INCEVT',PERnams));
        if any(ssm.INCEVT)==0
            try
                ssm.CEVT=CEVT;  % CEVT must be defined in mf_adapt 
            catch ME
               error('If any IN%s in the PER worksheet = 0 then you must define %s in mf_adapt !','CEVT');
            end
            for iComp=1:ssm.NCOMP
              cevt=PERvals(:,strmatchi(strmatchi('CEVT_%d',iComp),PERnams));
              for iPer=1:ssm.NPER
                  if ssm.INCEVT(iPer)>0, ssm.CEVT{iPer,iComp}=cevt(iPer); end
              end
            end
        else
            ssm.CEVT=cell(ssm.NPER,ssm.NCOMP);
            for iComp=1:ssm.NCOMP
                cevt=PERvals(:,strmatchi(sprintf('CEVT_%d',iComp),PERnams));
                for iPer=1:ssm.NPER
                    ssm.CEVT{iPer,iComp}=cevt(iPer);
                end
            end
        end
    end
        
    % point sources, are only specified in mf_adapt
    
    if exist('PNTSRC','var')
        if iscell(PNTSRC)
            ssm.PNTSRC =  PNTSRC;
        else
            ssm.PNTSRC = {PNTSRC};
        end
    end
        
    % count total number of sources and sinks in any one period included in
    % the flow model
    ssm.MXSS = MXSS + sum(IBOUND(:)<0 | ICBUND(:)<0);

    %% Generate point sources directly from WEL and MNW1

    %% Well package on?
    if namSCEN('WEL',nam.PCKG,nam.SCEN)

        % find wellObj's in the workspace
        I = strmatchi('wellObj',{mf_variables.class}); % <-- uses workspace variables

        if I
            % Initialize ssm.well as empty if it does not yet exist
            if ~isfield(ssm,'well'), ssm.well=[]; end

            % contatenate all wells in the workspace into ssm.well
            for i=1:numel(I)
                eval(['ssm.well = [ssm.well; ', mf_variables(I(i)).name, ']']); % <--- workspace variables
            end
            % remove wells without Q filled in.
            for iw = numel(ssm.well):-1:1
                if isempty(ssm.well(iw).Q) || all(isnan(ssm.well(iw).Q))
                    ssm.well(iw)=[];
                end
            end
        end
    end
        
    if namSCEN('MNW1',nam.PCKG,nam.SCEN)
  
        % find MNW1obj's in the workspace
        I = strmatchi('MNW1Obj',{mf_variables.class}); %<--- workspace variables
        if any(I)
            % initialize ssm.NMN1 as empty if it does not yet exist
            if ~isfield(ssm,'MNW1'), ssm.MNW1=[]; end

            % concatenate the MNW1's in the workspase into ssm.MNW1
            for i=1:numel(I)
                eval(['ssm.MNW1 =[ssm.MNW1; ',mf_variables(I(i)).name,']']); % <--- workshpace variables
            end
            % Only accept multinode wells with C filled in
            for iw =numel(ssm.MNW1):-1:1
                if isempty(ssm.MNW1(iw).C) || all(isnan(ssm.MNW1(iw).C(:)))
                    ssm.MNW1(iw)=[];
                end
            end
            % Verify that wll.C, well.Q and well.Dt are filled in
        end
    end
    
    if namSCEN('MNW2',nam.PCKG,nam.SCEN)
          
        % find MNW2obj's in the workspace
        I = strmatchi('MNW2Obj',{mf_variables.class});  % < --- workspace variables
        if I
            % initialize ssm.NMN2 as empty if it does not yet exist
            if ~isfield(ssm,'MNW2'), ssm.MNW2=[]; end

            % concatenate the MNW2's in the workspase into ssm.MNW1
            for i=1:numel(I)
                eval(['ssm.MNW2 =[ssm.MNW2; ',mf_variables(I(i)).name,']']); % <--- workspace variables
            end
            
            % Verify that wll.C, well.Q and well.Dt are filled in
            % Only accept multinode wells with C filled in
            for iw =numel(ssm.MNW2):-1:1
                if isempty(ssm.MNW2(iw).C) || all(isnan(ssm.MNW2(iw).C))
                    ssm.MNW2(iw)=[];
                end
            end
        end
    end
    
    
    %% add pointObj, lineObj and area2Obj
    % we already found the pointObjects, lineObjects and areaObjects when
    % preparing BCN above

    if exist('pointObjects','var') && ~isempty(pointObjects)
        ssm.point =  pointObjects(ismember({pointObjects.type},activeStresses));
    end
    if exist('lineObjects','var') && ~isempty(lineObjects)
        ssm.line  = lineObjects(ismember({lineObjects.type},activeStresses));
    end
    if exist('areaObjects','var') && ~isempty(areaObjects)
        ssm.area = areaObjects(ismember({areaObjects.type},activeStresses));
    end
    
    writeSSM(basename,ssm);
end

%% ===== THE GCG-file (Generalize conjugate gradient solver package)
gcg.SCEN=namSCEN('GCG',nam.PCKG,nam.SCEN);

if gcg.SCEN
    fprintf('Generating Generalized Conjugate Gradient Solver Process struct\n');
    gcg.unit=nam.UNIT(strmatchi('GCG',nam.PCKG));
    gcg.ext =nam.EXT {strmatchi('GCG',nam.PCKG)};
    %3
    gcg.MXITER=MT3Dvals(strmatchi('MXITER',MT3Dnams),gcg.SCEN);
    gcg.ITER1 =MT3Dvals(strmatchi('ITER1' ,MT3Dnams),gcg.SCEN);
    gcg.ISOLVE=MT3Dvals(strmatchi('ISOLVE',MT3Dnams),gcg.SCEN);
    gcg.NCRS  =MT3Dvals(strmatchi('NCRS'  ,MT3Dnams),gcg.SCEN);
    gcg.ACCL  =MT3Dvals(strmatchi('ACCL'  ,MT3Dnams),gcg.SCEN);
    gcg.CCLOSE=MT3Dvals(strmatchi('CCLOSE',MT3Dnams),gcg.SCEN);
    gcg.IPRGCG=MT3Dvals(strmatchi('IPRGCG',MT3Dnams),gcg.SCEN);

    writeGCG(basename,gcg);
end


%% ===== THE rct-file (chemical reaction package) =======
rct.SCEN=namSCEN('RCT',nam.PCKG,nam.SCEN);

% The general idea is that if a parameter is specified in the workspace of
% Matlab, then that parameter is used instead of the values specified in
% the workbook worksheet LAY. Note that we can only specify one value per
% layer in the workbook. To specify a value per cell, specify the parameter
% in Matlab's workspace as a 3D array or as a cell array in which each cell
% is a 3D array pertaining to the species in question.

if rct.SCEN
    rct.FREE =FREE;
    rct.NCOMP=btn.NCOMP;
    rct.MCOMP=btn.MCOMP;  % Don't know if this is needed in rct ## check
    rct.NCOL=GRID.Nx;
    rct.NROW=GRID.Ny;
    rct.NLAY=GRID.Nlay;
        
    fprintf('Generating Chemical Reaction Process struct\n');
    rct.unit=nam.UNIT(strmatchi('RCT',nam.PCKG));
    rct.ext =nam.EXT {strmatchi('RCT',nam.PCKG)};
    %E1 ISOTHM IREACT IRCTOP IGETSC
    %  ISOTHM sorption type
    %    0=no sorption
    %    1=linear sorption
    %    2=freundlich isotherm
    %    3=langmuir
    %    4=first order kinetic sorption (non-equilibrium)
    %    5=dual domain mass transfer (without sorption)
    %    6=dual domain mass transfer (with sorption)
    rct.ISOTHM=MT3Dvals(strmatchi('ISOTHM',MT3Dnams),rct.SCEN);

    % IREACT type of kinetic rate reaction flag
    %    0=no kinetic rate reaction
    %    1=first-order kinetic rate reaction
    %      chemical reaction simulations need add-on package
    rct.IREACT=MT3Dvals(strmatchi('IREACT',MT3Dnams),rct.SCEN);
  
    % IRCTOP reaction variable entry method flag
    %  >=2 all variables are entered as 3D array (using RARRAY)
    %  < 2 all variables are entered as 1D with one value per layer
    rct.IRCTOP=MT3Dvals(strmatchi('IRCTOP',MT3Dnams),rct.SCEN);
    if rct.IRCTOP<2
        error('IRCTOP = %d, this is deprecated and will produce no output! Set IRCTOP=2 in MT3D worksheet!',...
            rct.IRCTOP);
    end
    
    % IGETSC initial conc for the adsorbed phase reading flag (ISOTHM=4,5,6)
     rct.IGETSC=MT3Dvals(strmatchi('IGETSC',MT3Dnams),rct.SCEN);
     
     rct.NLAY=size(LAYparvals,1);
                    
     % Get columns form worksheet, note that these may be multiple
     % columns because we use more than one species
     if exist('RHOB','var')
         rct.RHOB = RHOB;
     else
         rct.RHOB = LAYparvals(:,strmatchi('RHOB'   ,LAYparnams));
     end
     
     if isvector(rct.RHOB) && GRID.Nlay>1, rct.RHOB=XS(rct.RHOB); end     
     if GRID.AXIAL, rct.RHOB = bsxfun(@times,GRID.TWOPIR,rct.RHOB); end

     if exist('PRSITY2','var')
         rct.PRSITY2 = PRSITY2;
     else
         rct.PRSITY2 = LAYparvals(:,strmatchi('PRSITY2',LAYparnams));
     end
     
     if isvector(rct.PRSITY2) && GRID.Nlay>1, rct.PRSITY2=XS(rct.PRSITY2); end         
     if GRID.AXIAL, rct.PRSITY2 = bsxfun(@times,GRID.TWOPIR,rct.PRSITY2); end

     if exist('SRCONC','var')
         if ~iscell(SRCONC)
             rct.SRCONC = {SRCONC};
         else
             rct.SRCONC = SRCONC;
         end
     else
         srconc = LAYparvals(:,strmatchi('srconc',LAYparnams));
         for i=1:rct.NCOMP
             rct.SRCONC{i}=XS(srconc(:,i));
         end
     end
     
     if exist('SP1','var')
         if isvector(SP1) && GRID.Nlay>1, SP1=XS(SP1); end
         if ~iscell(SP1)
             rct.SP1 = {SP1};
         else
             rct.SP1 = SP1;
         end
     else
         sp1 = LAYparvals(:,strmatchi('sp1',LAYparnams));
         for i=1:rct.NCOMP
             rct.SP1{i} = XS(sp1(:,i));
         end
     end
     
     if exist('SP2','var')
         if isvector(SP2) && GRID.Nlay>1, SP2 = XS(SP2); end
         if ~iscell(SP2)
             rct.SP2={SP2};
         else
             rct.SP2 = SP2;
         end
     else
         sp2 = LAYparvals(:,strmatchi('sp2',LAYparnams));         
         for i=1:rct.NCOMP
             rct.SP2{i}=XS(sp2(:,i));
         end
     end
     
     if exist('RC1','var')
         if isvector(RC1) && GRID.Nlay>1, RC1=XS(RC1); end
         if ~iscell(RC1)
             rct.RC1={RC1};
         else
             rct.RC1 = RC1;
         end
     else
         rc1 = LAYparvals(:,strmatchi('rc1',LAYparnams));         
         for i=1:rct.NCOMP
             rct.RC1{i}=XS(rc1(:,i));
         end
     end
     
     if exist('RC2','var')
         if isvector(RC2) && GRID.Nlay>1, RC2=XS(RC2); end
         if ~iscell(RC2)
             rct.RC2 = {RC2};
         else
             rct.RC2 = RC2;
         end
     else
         rc2 = LAYparvals(:,strmatchi('rc2',LAYparnams));
         for i=1:rct.NCOMP
             rct.RC2{i}=XS(rc2(:,i));
         end
     end
        
    writeRCT(basename,rct);
end

% ===== END of MT3DSM input ==============================

%% ======================== SWI input ============================
swi.SCEN=namSCEN('SWI',nam.PCKG,nam.SCEN,'exact');

%% Note that SWI computes freshwater heads in the top of the aquifers !

if swi.SCEN
    fprintf('Generating Salt Water Intrusion struct\n');
    swi.unit=nam.UNIT(strmatchi('SWI',nam.PCKG));
    swi.ext =nam.EXT {strmatchi('SWI',nam.PCKG)};
    swi.FREE     = FREE;
    swi.GRID     = GRID;    
    swi.ISWIZT   = MFLOWparvals(strmatchi('ISWIZT'  , MFLOWparnams),swi.SCEN);
    swi.ISTRAT   = MFLOWparvals(strmatchi('ISTRAT'  , MFLOWparnams),swi.SCEN);
    swi.NPRN     = MFLOWparvals(strmatchi('NPRN'    , MFLOWparnams),swi.SCEN);
    swi.TOESLOPE = MFLOWparvals(strmatchi('TOESLOPE', MFLOWparnams),swi.SCEN);
    swi.TIPSLOPE = MFLOWparvals(strmatchi('TIPSLOPE', MFLOWparnams),swi.SCEN);
    swi.ZETAMIN  = MFLOWparvals(strmatchi('ZETAMIN' , MFLOWparnams),swi.SCEN);
    swi.DELZETA  = MFLOWparvals(strmatchi('DELZETA' , MFLOWparnams),swi.SCEN);
    swi.NU       = MFLOWparvals(strmatchi('NU'      , MFLOWparnams),       :);
    try
        swi.ZETA = ZETA;
    catch ME
        error('mfLab:mf_setup:ZETA_missing_in_workspace',...
            ['You must define ZETA in the workspace if you use SWI !\n',...
             'Note that ZETA must either have the structure like that generated by readBud()\n',...
             'or be an Ny*Nx*NSRF 3D array in which each layer represents a ZETAPLANE in 3D space\n',...
             'in the latter case, no inversions are possible.\n']);
    end
    
    % NSRF = number of (concentration zones - 1)
    if isstruct(ZETA)
        if     isfield(ZETA,'values'), fieldNm = 'values';
        elseif isfield(ZETA,'term'),   fieldNm = 'term';
        else
            error('fieldNm of ZETA must be values or term');
        end
        swi.NSRF = numel(ZETA.(fieldNm)); % first item of arraystruct ZETA is used
    elseif isnumeric(ZETA)
        swi.NSRF = size(ZETA,3);
    else
        error('mfLab:mf_setup:SWI_wrong_structure',...
            ['ZETA for SWI must either be a struct(array) like that generated by readbud()\n',...
             'with field term{NSRF}(Ny,Nx,Nz) defining the ZETAPLANES,\n',...
             'or a Ny*Nx*NSRF array where each layer ipln defines the ZETAPLANE in 3D space.\n']);
    end
    
    try
        swi.SSZ = SSZ;     % SWI uses SSZ in place of PEFF
    catch ME
        try
            swi.SSZ = PEFF;    % SWI uses SSZ in place of PEFF
        catch ME
            error('mfLab:mf_setup:SWI_missing_PEFF',...
                'You must define SSZ or PEFF (=porosity) in the workspace \nwhen using SWI !');
        end
    end
    
    try
        swi.ISOURCE = ISOURCE;
        
        % If ISOURCE> 0 ? Sources and sinks are of the same type as water in zone ISOURCE.
        %      If such a zone is not present in the cell, sources and sinks are placed
        %      in zone at the top of the aquifer.
        % If ISOURCE= 0 ? Sources and sinks are of the same type of water as
        %       at the top of the aquifer.
        % If ISOURCE< 0 ? Sources are of the same type as water in zone ISOURCE.
        % Sinks are of the same type of water as at the top of the aquifer.
        % This option is useful for the modeling of the ocean bottom where
        % infiltrating water is salt, yet exfiltrating water is of the same
        % type as the water at the top of the aquifer.
    catch ME
        error('mfLab:ms_setup:ISOURCE_missing_in_workspace',...
            ['ISOURCE missing in workspace. ISOURCE(Ny,Nx,Nlay) required by SWI package\n',...
            ' ISOURCE>0 - Sources and sinks same type as water in zone ISOURCE\n',...
            ' ISOURCE=0 - Sources and sinks same type as water in top of aquifer\n',...
            ' ISOURCE<0 - Sources same type as water in zone ISOURCE, but\n',...
            '             Sinks same type as water in top of aquifer.\n',...
            '             Useful for modeling sea bottom on top of aquifer.\n']);
    end
   
    swi.NZONES = swi.NSRF+1;

    if swi.ISTRAT==1,
        swi.NU = swi.NU(1:swi.NSRF+1);
    else
        swi.NU = swi.NU(1:swi.NSRF+2);
    end
    
    writeSWI(basename,swi);
end


%% ======================== SWI2 input ============================
swi2.SCEN=namSCEN('SWI2',nam.PCKG,nam.SCEN);

%% Note that SWI2 computes freshwater heads in the top of the aquifers !

if swi2.SCEN
        
    [SWI2parnams,SWI2vals]=getExcelData(XLSF,'SWI2','Vertical');

    fprintf('Generating Salt Water Intrusion struct\n');
    swi2.unit=nam.UNIT(strmatchi('SWI2',nam.PCKG,'exact'));
    swi2.ext =nam.EXT {strmatchi('SWI2',nam.PCKG,'exact')};
    swi2.FREE     = FREE;
    swi2.GRID     = GRID;    
    
    %% first get ZETA and count number of surfaces NSURF
    try
        swi2.ZETA = ZETA;
    catch ME
        error('mfLab:mf_setup:ZETA_missing_in_workspace',...
            ['You must define ZETA in the workspace if you use SWI !\n',...
             'Note that ZETA must either have the structure like that generated by readBud()\n',...
             'or be an Ny*Nx*NSRF 3D array in which each layer represents a ZETAPLANE in 3D space\n',...
             'in the latter case, no inversions are possible.\n']);
    end
    
    if isstruct(ZETA)
        if     isfield(ZETA,'values'), fieldNm='values';
        elseif isfield(ZETA,'term')  , fieldNm='term';
        else
            error('if ZETA is a struct, fieldNm must be values or term');
        end
        swi2.NSRF = numel(ZETA.(fieldNm)); % first item of arraystruct ZETA is used

    elseif isnumeric(ZETA)
        swi2.NSRF = size(ZETA,3);
    
    else
        error('mfLab:mf_setup:SWI_wrong_structure',...
            ['ZETA for SWI must either be a struct(array) like that generated by readbud()\n',...
             'with field term{NSRF}(Ny,Nx,Nz) defining the ZETAPLANES,\n',...
             'or a Ny*Nx*NSRF array where each layer ipln defines the ZETAPLNAE in 3D space.\n']);
    end

    %% Then get the other variables
    
    swi2.ISTRAT   = SWI2vals(strmatchi('ISTRAT'  , SWI2parnams),swi2.SCEN);
    swi2.NOBS     = 0;
    swi2.ISWIZT   = SWI2vals(strmatchi('ISWIZT'  , SWI2parnams),swi2.SCEN);
    swi2.ISWIBD   = SWI2vals(strmatchi('ISWIBD'  , SWI2parnams),swi2.SCEN);
    swi2.ISWIOBS  = 0;
    swi2.ADAPTIVE = SWI2vals(strmatchi('ADAPTIVE' , SWI2parnams),swi2.SCEN);
    
    swi2.NSOLVER  = SWI2vals(strmatchi('NSOLVER' , SWI2parnams),swi2.SCEN);
    swi2.IPRSOL   = SWI2vals(strmatchi('IPRSOL'  , SWI2parnams),swi2.SCEN);
    swi2.MUTSOL   = SWI2vals(strmatchi('MUTSOL'  , SWI2parnams),swi2.SCEN);

    if swi2.NSOLVER ==2
        swi2.MXITER = SWI2vals(strmatchi('MXITER'  , SWI2parnams),swi2.SCEN);
        swi2.ITER1  = SWI2vals(strmatchi('ITER1'   , SWI2parnams),swi2.SCEN);        
        swi2.NPCOND = SWI2vals(strmatchi('NPCOND'  , SWI2parnams),swi2.SCEN);
        swi2.ZCLOSE = SWI2vals(strmatchi('ZCLOSE'  , SWI2parnams),swi2.SCEN);
        swi2.RCLOSE = SWI2vals(strmatchi('RCLOSE' , SWI2parnams),swi2.SCEN);
        swi2.RELAX  = SWI2vals(strmatchi('RELAX'  , SWI2parnams),swi2.SCEN);
        swi2.NBPOL  = SWI2vals(strmatchi('NBPOL'  , SWI2parnams),swi2.SCEN);
        swi2.DAMP   = SWI2vals(strmatchi('DAMP'   , SWI2parnams),swi2.SCEN);
        swi2.DAMPT  = SWI2vals(strmatchi('DAMPT'  , SWI2parnams),swi2.SCEN);
    end
    
    swi2.TOESLOPE = SWI2vals(strmatchi('TOESLOPE', SWI2parnams),swi2.SCEN);
    swi2.TIPSLOPE = SWI2vals(strmatchi('TIPSLOPE', SWI2parnams),swi2.SCEN);
    swi2.ALPHA    = SWI2vals(strmatchi('ALPHA'   , SWI2parnams),swi2.SCEN);
    swi2.BETA     = SWI2vals(strmatchi('BETA'    , SWI2parnams),swi2.SCEN);
    
    swi2.NADPTMX = SWI2vals(strmatchi('NADPTMX', SWI2parnams),swi2.SCEN);
    swi2.NADPTMN = SWI2vals(strmatchi('NADPTMN', SWI2parnams),swi2.SCEN);
    swi2.ADPTFCT = SWI2vals(strmatchi('ADPTFCT' , SWI2parnams),swi2.SCEN);

    if swi2.ISTRAT == 0
        swi2.NU       = SWI2vals(strmatchi('NU'      , SWI2parnams), 1:swi2.NSRF+2);
    else
        swi2.NU       = SWI2vals(strmatchi('NU'      , SWI2parnams), 1:swi2.NSRF+1);
    end
    
    try
        swi2.SSZ = SSZ;     % SWI uses SSZ in place of PEFF
    catch ME
        try
            swi2.SSZ = PEFF;    % SWI uses SSZ in place of PEFF
        catch ME
            error('mfLab:mf_setup:SWI_missing_PEFF',...
                'You must define SSZ or PEFF (=porosity) in the workspace \nwhen using SWI !');
        end
    end
    
    try
        swi2.ISOURCE = ISOURCE;        
        % If ISOURCE> 0 ? Sources and sinks are of the same type as water in zone ISOURCE.
        %      If such a zone is not present in the cell, sources and sinks are placed
        %      in zone at the top of the aquifer.
        % If ISOURCE= 0 ? Sources and sinks are of the same type of water as
        %       at the top of the aquifer.
        % If ISOURCE< 0 ? Sources are of the same type as water in zone ISOURCE.
        % Sinks are of the same type of water as at the top of the aquifer.
        % This option is useful for the modeling of the ocean bottom where
        % infiltrating water is salt, yet exfiltrating water is of the same
        % type as the water at the top of the aquifer.
    catch ME
        error('mfLab:ms_setup:ISOURCE_missing_in_workspace',...
            ['ISOURCE missing in workspace. ISOURCE(Ny,Nx,Nlay) required by SWI package\n',...
            ' ISOURCE>0 - Sources and sinks same type as water in zone ISOURCE\n',...
            ' ISOURCE=0 - Sources and sinks same type as water in top of aquifer\n',...
            ' ISOURCE<0 - Sources same type as water in zone ISOURCE, but\n',...
            '             Sinks same type as water in top of aquifer.\n',...
            '             Useful for modeling sea bottom on top of aquifer.\n']);
    end
    
    writeSWI2(basename,swi2);
end

%% ============ End of input ============
fprintf('\nM O D E L   S U M M A R Y\n');
fprintf('Packages that are on:\n');
fprintfs(1,'   %s',nam.PCKG); fprintf('\n');
fprintf('Model limits and size:\n');
fprintf('xLim = [%g %g], Lx = %g\n',gr.xGr([1 end]),diff(gr.xGr([1 end])));
fprintf('yLim = [%g %g], Ly = %g\n',gr.yGr([end 1]),diff(gr.yGr([1 end])));
fprintf('zMax = %g, zMin = %g\n',gr.zGr([1 end]));
fprintf('Model grid:\n');
fprintf('GRID = [%d %d %d] MINDZ=%g, LAYCBD =',GRID.size,GRID.MINDZ);
if all(GRID.LAYCBD==0)
    fprintf(' LAYCBD: all zero\n');
else
    fprintf(' LAYCBD:');
    fprintf(' %d',GRID.LAYCBD); fprintf('\n');
end
fprintf('Model type:\n');
fprintf('AXIAL = %d\n',GRID.AXIAL);

fclose('all');  % just to make sure all files are closed
fprintf('.... finished generation of input files, all files closed! Time= %s!\n\n',datestr(now));

toc;

%% ========= RUN model(s): Seawat MT3DMS or mf2k =====================

% Execute model code in background?
if BACKGROUND, BACKGROUND=' &'; else BACKGROUND=''; end
   
if ~isempty(GREP), GREP=['| grep "' GREP '"']; end

if ~ismac,   GREP='';  end

quote='"';

tic;

%% We run all selected models, MODFLOW is automatically set if MT3DMS is chosen
%  otherwise we run models in the order of the packages in the NAM sheet
%  so MODPATH must come after MODDLOW
NModel=numel(nam.MODEL);

for iModel=1:NModel
        
    if iModel<NModel, RUNBACKGROUND=''; end
    
    [MP,MF,FE]=fileparts(nam.mdlpath{iModel});
    
    if ~exist(nam.mdlpath{iModel},'file'),
        error(['Sorry, mfLab can''t find find executable <<%s>>.\n',...
            'Most probably this is because the model you have chosen is only available under Windows.\n',...
            'Check your NAM worksheet.'],nam.mdlpath{iModel});
    end

    switch lower(nam.MODEL{iModel})
        case 'whatever', %'MFSWI'
            [status,result]=system([quote,nam.mdlpath{iModel},quote,' ',MF,'.nam  ']); % runBckGrd does not work with SWI             fprintf('%3i: %s status = %d, result=%s\n',i,nam.MODEL{iModel},status,result);
            fprintf('%3i: %s status = %d, result=%s\n',i,nam.MODEL{iModel},status,result);
        case 'modpath'
            % Modpath6.0 expects binary files without record lenght info,
            % so-called streamed files (ACCESS/'STREAM'/ not ACCESS/'SEQUENTIAL'
            % cleanBinary removes this record info if present.
            cleanBinary([basename '.HDS']);
            cleanBinary([basename '.BGT']);
            call = [quote,nam.mdlpath{iModel},quote,' ',pth.simulationFile,' ',GREP,' ',BACKGROUND];
            fprintf('Trying to run \n%s\n',call);
            system(call);
        otherwise % any other model
            call = [quote,nam.mdlpath{iModel},quote,' ',MF,'.nam',' ',GREP,' ',BACKGROUND];
            fprintf('Trying to run \n%s\n',call);
            system(call);
            if iModel<NModel,
                wait4model(MF);
            end
    end
end
    

if ~isempty(BACKGROUND)
    cprintf('blue','\nWARNING:\n');
    cprintf('blue','Model code is now running in the background!\n');
    cprintf('blue','Running in the background allows your to continue using Matlab wihtout waiting for the results\n\n');
    cprintf('blue','You may change the variable BACKGROUND in mf_adapt/mf_build to 0 to prevent that!\n');
    cprintf('blue','While running in the background you may run sytem commands like\n');
    cprintf('blue','system(''tail -100 *.LST | grep SUMMARY''); to see the tail of the LST file (UNIX, MAC)\n');
    cprintf('blue','To see the running processes in the task window of MS Windows press CTRL-ALT-DEL\n');
    cprintf('blue','You may load the list (.LST) file in a good editor and watch it grow during the run.\n');
    cprintf('blue','Or just continue working until the fan becomes quiet, indicating the run has finished.\n');
    fprintf('\n');
    cprintf('blue','Make sure you wait until the model finished before running mf_adapt !!!\n\n');
end

%%  give the user its grid back

fprintf('Time to run model is %g seconds\n',toc);

%% now run mf_analyze

%% Set AFTERMFSETUP to run this file at the end of mf_setup, used for calibration

if exist('AFTERMFSETUP','var') && ~strcmp(runBckGrd,' &')
    eval(AFTERMFSETUP);
end
