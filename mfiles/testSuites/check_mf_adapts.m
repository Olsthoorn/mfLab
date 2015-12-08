function check_mf_adapts(actualPath)
%CHECK_MF_ADAPTS runs mf_adapt.m and mf_build.m in current and all subdirectories
%
% Example:
%    check_mf_adapts();
%
% What this function does:
%    1) visit each subdiretory recursively
%    2) test if mf_adapt exists
%    3) run it, verify there is no failure
%    4) then retract and do the next
%
% TO 130401

global homeDir checkAgain

testInterval = 10; % days

if nargin==0
    homeDir     = pwd;
    checkAgain  = now+5;
    check_mf_adapts(homeDir);
    return;
end

%fprintf('running: %s\n',actualPath);

d = dir(actualPath);

d=d(~ismember(1:numel(d),strmatchi('.',{d.name})));

for id = find([d.isdir])
    fprintf('>');
    check_mf_adapts(fullfile(actualPath,d(id).name));
end

%% Alternative names for mf_adapt
names = {'mf_adapt.m', 'mf_build.m'};

for i = 1:numel(names)
    found = strmatchi(names{i},{d.name},'exact');
    if found
        [~,name] = fileparts(names{i});
  %      fprintf('%s found\n',name);        
        cd(actualPath) %#ok
        try load('testLog')
            mustTest = now>testTime+testInterval; %#ok
        catch ME
            mustTest = true;
            fprintf('%s\n',ME.message');
        end
    
        if mustTest
            save homeDir homeDir
            fprintf('running  %s\n',name);
            evalin('base',name);
            mf_cleandir
            close('all');
            testTime =now; %#ok from testLog
            save('testLog','testTime');
            load  homeDir
            delete homeDir.mat
        else
            fprintf('skipped %s\n',name);
        end
        break;
    end
end

fprintf('<\n');

cd(homeDir); %#ok
