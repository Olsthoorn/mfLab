function addcopyright(varargin)
%ADDCOPYRIGHT adds or updates mfLab copyright notice
%
% USAGE:
%   addcopyright(files) --- adds or updates mfLab copyright notice to mfiles.
%   addcopyright(dir)   --- adds or updates mfLab copyright notice to mfiles in dir
%   addcopyright()      --- adds or updates mfLab copyright notice to current dir
%
% TO 130427

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

cpr= sprintf('\n%% Copyright %d-%s Theo Olsthoorn, TU-Delft and Waternet, without any warranty\n%s',...
              2009,datestr(now,'yyyy'),'% under free software foundation GNU license version 3 or later');

mfileDirs = {'mfiles','write','read','etc','analytic','visualization','gridcoords','fdm','NHI','testSuites'};

%% only mfile directories of mfLab
[~,subDir] = fileparts(pwd);
if ~strmatchi(subDir,mfileDirs)
    error('%s: for this function move to one of the directories under mfLab/mfiles first',mfilename);
end

if isempty(varargin)
    d = dir('*.m');
else
    [~,~,e] = fileparts(varargin{1});
    if strcmp('.m',e)
        d = dir(varargin{1});
    else
        d = dir([varargin{1} '/*.m']);
    end
end

for i=1:length(d)

    %% do only m files
    if ~(d(i).name(end-1) == '.' && d(i).name(end)=='m')
        continue;
    end
        
    %% don't add copyright to certain m-files
    if strcmpi(d(i).name,{'mf_adapt.m','mf_build.m','mf_analyze.m'})
        continue;
    end
        
    fid=Fopen(d(i).name,'r');

    %% hoist entire file into string
    s = fscanf(fid,'%c',[1,Inf]);

    fclose(fid);
    
    %% find copyright notice specific of mfLab
    I = regexpi(s,'^%[ ]*Copyright[ ]*2009| Theo Olsthoorn');

    switch numel(I)
        case 2  % it's a hit, both copyright and Theo Olsthoorn found
            fprintf('Updating copyright notice of file %s\n',d(i).name);
            s = [s(1:I(1)) sprintf(' Copyright 2009-%s',datestr(now,'yyyy')) s(I(2):end)];
        case 0 % no copyright at all present, so add one
            fprintf('Adding copyright notice to file %s\n',d(i).name);

            % get lines starting with %
            J = regexp(s,'\n%'); if isempty(J), J=1; end

            % get lines ~starting with %
            K = regexp(s,'\n[^%]');

            % get first non-comment line after first comment line
            K = K(K>J(1));
            K = K(1);

            % merge copyright notice
            s = [s(1:K) cpr s(K:end)];
        otherwise
            % do nothing
            fprintf('File %s already has copyright notice to file %s\n',d(i).name);
            continue;
    end

    % saveguard the old mfile first

    movefile(  d(i).name, [d(i).name '_old']);
    
    % write s to new file
    fid = Fopen(d(i).name    ,'w');

    fprintf(fid,'%c',s);

    fclose(fid);

    fprintf('File %s done\n',d(i).name);
end
    
function fid=Fopen(fname,mode)
% fid=Fopen(fname,mode)
% fopen with error control

fid=fopen(fname,mode);
if fid<0
    error('Can''t open file <<%s>> for mode <<%s>>',fname,mode);
end
