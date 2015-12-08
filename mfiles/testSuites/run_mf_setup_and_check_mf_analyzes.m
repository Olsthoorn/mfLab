function run_mf_setup_and_check_mf_analyzes(varargin)
%RUN_MF_SETUP_AND_CHECK_MF_ANALYZES recursively check all examples
%
% Example:
%    run_mf_setup_and_check_mf_analyzes()
%
% What it does
%     run_mf_setup_and_check_mf_analyzes(actualPath,options)
%     recursively run mf_setup in every directory where an mf_adapt or an
%     mf_build is found:
%     In steps:
%       1) visit each subdiretory recursively
%       2) test if mf_adapt exists
%       3) run mf_setup
%       4) verify there is no failure
%       5) the run mf_analyze
%       6) if success,
%            then leave message with time step
%       7) retract and do the next
%
% Example:
%       run mf_setup_and_check_mf_analyzes();
%       run_mf_setup_and_check_mf_analyzes(actualPath)
%       run_mf_setup_and_check_mf_analyzes(actualPath,options)
%       run_mf_setup_and_check_mf_analyzes(options)
%       run_mf_setup_and_check_mf_analyzes('reCheck',now)
%
% options is a property,value pair. Currently the following combinations are understood:
%     'newDate',datenum
%          in all directories with testLog2.mat, replaces the checkdate
%          with a new value. This allows postponing the checkdate to the
%          future or to the past, whatever is required. Resetting
%          everything for testing the next time, can be done by setting the
%          newData to now().
%
% TO 120424

global homeDir checkAgain

% homeDir is the top of the tree to be searched and where to return to
% after the search completed sucessfully
% checkAgain is a matlab dateNum that will be stored in testLog2.mat in
% every directory with mf_adapt/mf_build and mf_analyze marking the moment
% after which the test will be run in case this directory is processed in
% the future.

[newDate,varargin] = getProp(varargin,{'newDate','reCheck','checkAgain'},[]);
if ~isempty(newDate)
    checkAgain = newDate;
else
    if isempty(checkAgain)
        checkAgain = now+30;
    end
end

if isempty(varargin)
    homeDir     = pwd;
    run_mf_setup_and_check_mf_analyzes(homeDir);
    return;
end

actualPath = varargin{1};

fprintf('%s\n',actualPath);

d = dir(actualPath);

d=d(~ismember(1:numel(d),strmatchi('.',{d.name})));

for id = find([d.isdir])
    fprintf('>');
    run_mf_setup_and_check_mf_analyzes(fullfile(actualPath,d(id).name));
end

%% Alternative names for mf_adapt
names = {'mf_adapt.m'; 'mf_build.m'};

for i = 1:numel(names)
    found = strmatchi(names{i},{d.name},'exact');
    if found
        [~,name] = fileparts(names{i});
        fprintf('%s found\n',name);        
        cd(actualPath) %#ok
        try load('testLog2')
            % testLog2 contains variables testTime and runTime
            if checkAgain<testTime
                testTime=checkAgain;
                save('testLog2','testTime','runTime');
                fprintf('TestTime reset to: %s in %s',testTime,pwd);
            end
            mustTest = now>testTime;
        catch ME
            mustTest = true;
            fprintf('%s\n',ME.message');
        end
    
        if mustTest
            save homeDir homeDir
            mf_cleandir
            save('adaptName','name');  % must save because of clearing in mf_setup
            save('checkTime','checkAgain')
            tic;
            mf_setup;
            load('adaptName'); delete('adaptName.mat');  % retrieves name
            if fileGrep([name '.m'],'BACKGROUND *= *1|BACKGROUND *= *true')
               wait4model(MF); % MF comes from mf_setup script (name of last model)
            end            
            checkHds(dir('*.HDS')); % test if succesfull
            if ~isempty(dir('mf_analyze.m'))
                mf_analyze;  % mf_analyze absent in tutorial
            end
            mf_cleandir
            close('all');
            load('checkTime'); delete('checkTime.mat')
            if checkAgain<now
                testTime = checkAgain + 30;
            else
                testTime = checkAgain; %#ok from testLog
            end
            runTime  = toc; %#ok
            save('testLog2','testTime','runTime');
            load  homeDir
            delete homeDir.mat
        end
        break;
    end
end

fprintf('<\n');

cd(homeDir); %#ok
