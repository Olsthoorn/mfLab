function look_for_use_of_function(funName,retest,actualPath)
%LOOK_FOR_USE_OF_FUNCTION checks in current and lower dirctories where string is used in m-files
%
% Example:
%    look_for_use_of_function(funName,retest);            % retest = false or true or empty, default = false
%    look_for_use_of_function('sinespace');               % where is sinespace used
%    look_for_use_of_function('sinespace',true)           % overwrite recheck date
%
% What it does:
%    1) visit current and all subdirectories recursively
%    2) look_for_use of string in local m-files
%    3) if found load mfile in editor and stop.
%    4) after editing press contiue on debug continue button above screen
%    5) return to home diretory when finished
%
% This function is useful to update calls of functions througout all directories
% as well as to look where a certain string is used in any mfile in the
% current and lower directories.
%
% TO 120323

global checkAgain homeDir

if nargin<1
    error('%s: must supply a string to look for in m files from given directory root',mfilename);
end

if nargin<2
    retest=false;
end

if nargin<3
    if ~islogical(retest)
        error('%s: use false or true for second argument',mfilename);
    end
    homeDir     = pwd;
    checkAgain  = now+5;
    
    if retest, RETEST='TRUE'; else RETEST='FALSE'; end
    
    fprintf(['\n\n%s:\nChecking this and all subdirectories recursively for mfiles\n',...
        'containing the string ''%s'' with retest option set to %s.\n'],...
        mfilename,funName,RETEST);
    
    look_for_use_of_function(funName,retest,homeDir);
    return;
end

fprintf('%s\n',actualPath);

d = dir(actualPath);

d=d(~ismember(1:numel(d),strmatchi('.',{d.name})));

for id = 1:numel(d)
    if d(id).isdir
        fprintf('>');
        look_for_use_of_function(funName,retest,fullfile(actualPath,d(id).name));
    else
        [~,fname,ext] = fileparts(d(id).name);
%        fprintf(1,'%s\n',fullfile(actualPath,fname,ext));
        if strcmp(ext,'.m')
            cd(actualPath) %#ok
            mfile = fullfile(actualPath,[fname,ext]);
            if ~isempty(dir('testLogFun.mat'))
                load('testLogFun')
                mustTest = retest || now>testTime;
            else
                mustTest = true;
            end
            if mustTest
                fid = fopen(mfile,'rt');
                if fid<0, error('%s: Can''t open file <<%s>>',mfile); end                
                if any(regexp(fscanf(fid,'%c',[1,Inf]),funName));
                    testTime = checkAgain;                    
                    save('testLogFun','testTime');
                    beep;
                    edit(mfile);
                    eval('keyboard');
                else
                    fclose(fid);
                end
            end
        else
            % skip file
        end
    end
end

fprintf('<\n');

cd(homeDir); %#ok
