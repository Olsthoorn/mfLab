function seeIfEverythingWorks(fp,thisDir)
% seeIfEveryThingWorks(); % Verifies that all example model work.
%
% From a given nodal directory looks hierarchially into all subfolders and
% when arrived at a leaf checks existance of
%   mf_adapt.m | mf_build.m
%   <<basename>>.xls
%   mf_analyze.m
%
%  where <<basename>> is scanned from the mf_adapt.m or mf_build.m.
%
% Wheather these file were found or not found is written to the log file,
% called mainlog-yyyy-mm-dd-hh-MM-ss.txt in the nodal start directory.
%
% If the for mfLab required files mf_adapt|mf_build, mf_analyze and
% <<basename>>.xls exist, it will run
%    mf_adapt or mf_build
%    mf_setup
%    mf_analyze
% it will capture the exit status of each of these procedures and and
% report that to the log file.
% It then continues, untill the last directory has been processed.
%
% If modflow, mt3dms or seawat files exist in the diretory, mf_setup is not
% run, only mf_adapt|mf_build and mf_analyze are run.
%
% if mf_setup was run, it will clearout the model files afterwards to
% prevent clutter of the hard disk.
%
% TO 121018

%% Initialize at startup with empy arguments
if nargin==0,
    fp = fopen(['mainlog-' datestr(now,'yyyy-mm-dd-hh-MM-ss') '.txt'],'w');
    thisDir = [pwd filesep];
end

%% get supbdirs excluding current and parent
d = subdirs(thisDir);
if numel(d)>0
    % for each subdir see if everything works (starging at the deepest)
    for iDir = numel(d):-1:1
        thatDir = [thisDir d(iDir).name '/'];
        seeIfEverythingWorks(fp,thatDir);
    end
else
    % if nosubdirs verify current one
        
    if  exist([thisDir 'mf_adapt.m'],'file')
        basename = getBasename([thisDir 'mf_adapt.m']);
        ada = 'mf_adapt.m';
    elseif exist([thisDir 'mf_build.m'],'file');
        basename = getBasename([thisDir 'mf_build.m']);
        ada = 'mf_build.m';
    else
        basename  = '';
        ada = 'NO-mf_adapt.m';
    end
 
    if ~isempty(basename) && exist([thisDir basename '.xls'],'file')
        xls = [basename '.xls'];
    else
        basename = 'NO-basename';
        xls = [basename '.xls'];
    end

    
    if exist([thisDir 'mf_analyze.m'],'file')
        ana = 'mf_analyze';
    else
        ana = 'NO_mf_analyze';
    end
    
    if ~strcmpi(ada,'NO-mf_adapt.m')   && ... 
       ~strcmpi(ada,'NO-mf-build.m')   && ...
       ~strcmpi(ana,'NO-mf_analyze.m') && ...
       ~strcmpi(xls,'NO-basename.xls')
   
   % running mf_adapt or mf_build
        fprintf('Running  %s\n',[thisDir ada]);
        
    % running mf_setup if necessary
        d = dir([thisDir '*.bat']);
        if isempty(d)
            fprintf('Running  mf_setup with %s\n',[thisDir xls]);
        else
            fprintf('Skipping mf_setup with %s\n',[thisDir xls]);
        end
        
    % running mf_analyze
        fprintf('Running  %s\n',[thisDir ana]);
        
    % cleaning up
        if isempty(d)
            mf_cleandir(thisDir);
        end
    end
    
    fprintf(fp,'%s    %-20s %-20s %-20s ''%s''\n',datestr(now),ada,ana,xls,thisDir);
end

if nargin==0,
    fclose(fp);
end

end

function d = subdirs(path)
    d=dir(path);
    d= d([d.isdir]);
    
    I = true(size(d));
    I(strmatchi('.',{d.name}))=false;
    
    d=d(I);
end

function basename = getBasename(fname)
    fp = fopen(fname,'r');
    
    if fp<1, error('%s: Can''t open file <<%s>>',mfilename,fname); end
        
    basename = fscanf(fp,'%c');
    basename = basename(regexp(basename,'basename*='):end);   
    basename = basename(regexp(basename,'='):regexp(basename,';'));
    basename = basename(2:end-1);
    basename(basename==' ')=[];
    basename(basename=='''')=[];

    fclose(fp);
end
