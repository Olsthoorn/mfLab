function mf_cleandir(thisDir)
%MF_CLEANDIR safely removes the files named in the nam files in this directory and some more
%
% USAGE:
%    mf_cleandir([dir]);  
%
% files with these extensions are also removed
% 'bat','mat','for','p00','err','UCN','CNF','MAS','BCF','OBS'
%
%   TO 091127 110514 121018
%
% Example:
%   mf_cleandir
%   mf_cleandir(dir);
%
%
% Copyright 2009-2011 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if nargin==0, thisDir=pwd; end

fprintf('Cleaning %s\n',thisDir);

%% delete all files thata are specified in the .nam files
if nargin==0,  thisDir =[pwd '/']; end
    
d=dir([thisDir '*.nam']);

for i=1:length(d)
    fid=fopen([thisDir d(i).name],'r');
    
    fgetl(fid); s=fgetl(fid);
    while s~=-1
        s     = s(find(s==' ',1,'last')+1:end);
        fName = [thisDir s];
        if exist(fName,'file')
            fprintf('Deleting %s\n',fName);
            delete(fName);
        end
        s=fgetl(fid);
    end
    fclose(fid);
end

%% some matfiles can always be removed
expr ='underneath|yourgridname|CINACT|HINACT|name';

d = dir([thisDir '*.mat']);
for id=1:length(d)
    if ~isempty(regexp(d(id).name,expr,'ONCE'))
        fprintf('Deleting %s\n',[thisDir d(id).name]);
        delete([thisDir d(id).name]);
    end
end

%% file with certain extension can always be removed (not in the nam files)
ext={'nam','bat','m~','for','p00','err','UCN','CNF','MAS','BCF','OBS','FTL','MNWa','MNWb','MNWc'};

for ie=1:length(ext)
    d=dir([thisDir '*.' ext{ie}]);
    for i=1:length(d)
        if strcmp(ext{ie},'mat')
            R=input(sprintf('\nSure to delete %s ??  [Y/N]',d(i).name),'s');
            if ~isempty(R) && lower(R(1))=='y',
                fprintf('Deleting %s\n',[thisDir d(i).name]);
                delete([thisDir d(i).name]);
            else
                fprintf('Skipping %s\n',[thisDir d(i).name]);
            end
        else
            fprintf('Deleting %s\n',[thisDir d(i).name]);
            delete([thisDir d(i).name]);
        end
    end
end

if exist('Underneath.mat','file'), delete('Underneath.mat'); end

delete('OUTput_*');  % MNW files

%%% Remove mdopath files
delete('mp6.*');
delete('MPATH6.*');

%% apply regexp if expression is given (this is TOO dangerous)
% better specify positively what can be removed.
% d = dir(thisDir);
% d = d(~[d.isdir]);
% for id=1:length(d)
%     if  isempty(regexp(d(id).name,'\.m$|\.xls?|\.png|\.jpg|\.pic|\.avi|\.mat|\.txt|\.inc','ONCE'))  % never remove these whatever the user wants
%         delete([thisDir d(id).name]);
%         fprintf('Deleting %s\n',[thisDir d(id).name]);
%     else
%         fprintf('Skipping %s\n',[thisDir d(id).name]);
%     end
% end
