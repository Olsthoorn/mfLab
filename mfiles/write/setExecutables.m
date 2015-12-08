%SETEXECUTABLES script that sets mfLab's executables in mflab/bin
%
% This scipt is put in mflab\mfiles\write together with the file mf_setup.m
% which calls this script.
% This script, setExecutables.m, is called as a script by
% mf_setup to set the paths to the executables (mf2k, mt3dms, seawat etc)
% As in this example you may build in switches in case you use mfLab in
% different environments like the PC and the MAC. As in this case the same
% directories have different names viewed from the MAC or the Windows
% operating system on the same computer. This is taken care of by the
% switch cause by the if statement. The switch makes use of the function
% ismac and ispc. These are avaialble from Matlab version 7. If not you may
% look at the function computer. If you only work on the PC you may delete
% the mac portion altogether.
% It may be convenient to exclude this script from the svn version control,
% so that it will not be updated inadvertently. In that case future updates of
% mfLab should be guaranteed to work immediately because they don't affect
% the specific file locations on you computer. 

% TO 091218
% TO 110513(mf2007 (CPF)) added, may replace mf2005 in future if proven it
% contains everything of mf2005 plus CFP.

% Path to the executables. I put all of them into mfLab/MODELS/bin
% for convenience, because I have several differntly compiled versions.
% It's your choice however. Anyway, set the parameters MODFLOW MT3D etc
% down below to their actual lcoations on your hard drive.
%
% If you don't like to copy your executables to another location, then copy
% a link to them into the mflab/bin directory.

% NOTICE the " " in the paths below to manage spaces in file names
% in system command used to launch the external executable later on

fprintf('Defining paths to your executables\n');

if ismac
    MODELS='/Users/Theo/GRWMODELS/mflab/bin/';  % location of my executables
%    if verLessThan('matlab', '8.0.1')
%        MF2000 =[MODELS,'mf2k-i386.mac'    ];  % location of MODFLOW executable
%        MF2005 =[MODELS,'mf2005-i386.mac'  ];  % location of MODFLOW 2005 executable
%    else
        MF2000 =[MODELS,'mf2k.mac'    ];  % location of MODFLOW executable
        MF2005 =[MODELS,'mf2005.mac'  ];  % location of MODFLOW 2005 executable
%    end
    MF2007 =[MODELS,'mf2005cfp.mac'  ];  % location of MODFLOW 2007 executable (CFP)
    MF2005NWT = [MODELS 'mf2005nwt.mac']; % location of MODFLOW-NWT excutable
    MT3DMS =[MODELS,'mt3dms5s.mac'];  % MT3DMS executable
    SEAWAT =[MODELS,'swt_v4.mac'  ];  % SEAWAT Executable
    MFSWI  =[MODELS,'mf2kswi.mac' ]; % mf2005 which knows SWI
    MODPATH=[MODELS,'mp6.mac' ]; % modpath
elseif ispc
    MODELS='Z:\mflab\bin\'; % location of my executables
    MF2000 =[MODELS,'mf2k.exe'    ];  % location of MODFLOW executable
    MF2005 =[MODELS 'mf2005.exe'  ]; % location of MF2005 executable
    MF2007 =[MODELS 'mf2005cfp.exe']; % location of MF2005 executable (CFP)
    MF2005NWT = [MODELS 'mf2005nwt.exe']; % locatin of MODFLOW-NWT executable
    MF2KASP=[MODELS 'mf2kasp.exe' ];   % Doherty's mf2k with improved for dry cells
    MT3DMS =[MODELS,'mt3dms5b.exe'];  % MT3DMS executable (use binary version on windows
          % so that it is comopatible with the stanard windows mf2k exacutable)
    SEAWAT =[MODELS,'swt_v4.exe'  ];  % SEAWAT Executable
    MFSWI  =[MODELS,'mf2kswi.exe' ];  % mf2k which knows SWI
    MODPATH=[MODELS,'mp6.exe'];         % modpath version 6
else
    error(['You are running <<%s>>, i.e. on either a Mac or PC, as yet only mac and pc are supported\n',...
           'Nevertheless, unix is expected to run on mac without change.\n',...
           'Try changin ismac to isunix in setup, rarray and m files readding unformatted\n',...
           'model output (i.e. readDat, readBud, readMT3D).\n',...
           'Type help computer for more information.'],computer);
end
