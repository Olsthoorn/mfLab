%MFSETPATH  --- set mflab paths and provide strings to be copied into a shortcut
% run this from the diretory mflab/mfiles
% to set the paths for the current session
%
% for a more permanent fixation of the path each time a matlab session is
% opened, you best way possibly is to put these instructions in a shortcut
% on the shortcut bar (not that from Matlbab 2013 the user inteface of
% Matlab has changed, which is also the case with the shortcut bar).
%
% To extend the search path with the directories that hold the mfLab
% mfiles, copy the lines below in a new shortcut. But make sure that the
% top line, which sets the path, matches the actual path to the mfLab/mfiles
% dirctory on your computer. You can get the full path easily by navigating
% to the directory mfLab/mfiles and then typien pwd (print working
% directory] the path name can then be copied into the shortcut.
%
% Windows uses use backward slashes \, unix and Mac use forward slashes /
%
% Example: 
% mfsetpath
%
% Note: call mfsetpath from the mfLab/mfiles directory type mfsetpath
%
% See Also: setexecutables
%
% setexecutables containes the names of the excutables that mfLab
% recognizes. Have a look and make sure the paths in it refer to the
% directory holding the executables, i.e mfLab/bin

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

pathstr = [fileparts(which('mfsetpath')) filesep];

addpath([pathstr 'analytic']);
addpath([pathstr 'etc']);
addpath([pathstr 'fdm']);
addpath([pathstr 'NHI']);
addpath([pathstr 'gridcoords']);
addpath([pathstr 'read']);
addpath([pathstr 'visualization']);
addpath([pathstr 'write']);

fprintf('The mfLab paths have been set for the current sesion.\n\n');

%% Instructions to make it easier by creating a shortcut

fprintf('These are strings you may want to put in a shortcut to easily set\n')
fprintf('your paths next time by a single press. See the Shortcuts toolbar\n');
fprintf('on the top of the Matlab screen\n\n');

mdirs={'analytic','etc','fdm','gridcoords','read','visualization','write'};

for i=1:length(mdirs)
    fprintf('addpath(''%s'');\n',[pathstr,mdirs{i}]);
end

fprintf('\nHOW TO CREATE A SHORTCUT ?\n');
fprintf('1) Click the right-hand button when the mouse is in the shortcut toolbar\n');
fprintf('   just below the icons on the top left of the Matlab window.\n');
fprintf('2) Select New Shortcut to open\n');
fprintf('3) Give it a name, like mfpaths\n');
fprintf('4) Copy the "addpath" lines above into the large textbox and close it\n');
fprintf('\nNext time when mfLab cannnot find its files, just press this shortcut\n');
fprintf('\n');


%% filling in the pathtoexecutables in mflab/mfiles/write/setExecutables.m

%% get path to executables relative to current directory
cd('..');

fprintf('\nChanged directory to\n%s\n',pwd);

pathtoexecutables=[pwd filesep 'bin' filesep];
fprintf('\npath to executables is\n%s\n',pathtoexecutables);

cd('mfiles');
fprintf('\nChanged directory to\n%s\n',pwd);

%% Read setPathtoexecutables.txt

fid=fopen([pathstr 'write' filesep 'setExecutables.m'],'r');

A=fscanf(fid,'%c',Inf);

fclose(fid);

%% Writing setPathtoexecutables.m replacing the string pathtoexecutables here

fprintf('\n');
fprintf('Writing setPathtoexecutables.m replacing the string\n');
fprintf('pathtoexecutables\n');
fprintf('herein with the actual path to the executables, which usually is\n');
fprintf('Drive:\\<<path to mflab>>\\bin\n');
fprintf('If it doesn''t work you have to do it by hand.\n');
fprintf('On PC change the path to the form like\n');
fprintf('''%s%sbin''\n',pwd,filesep);
fprintf('instead of something like\n');
fprintf('''vmware-host\\Shared Folders\\tolsthoorn On My Mac\\GRWMODELS\\mflab\\bin''\n');
fprintf('otherwise windows is unaware of file changes during a session\n');
fprintf('\n');

fid=fopen([pathstr filesep 'write' filesep 'setExecutables.m'],'w');

n=strfind(A,'MODELS=')+length('MODELS=')-1;

if ismac, n=n(1); else n=n(2); end

m=strfind(A(n:end),''''); m=m(2);

fprintf(fid,'%s',A(1:n));

fprintf(fid,'%c%s%c','''',pathtoexecutables,'''');


fprintf(fid,'%s',A(n+m:end));

fclose(fid);
