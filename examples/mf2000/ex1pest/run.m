% This is the mfile used as argument when invoking matlab from the command
% line (by PEST or yourself) using the DOS command:
% matlab -nosplash -nodesktop -r run.n -logfile run.log
% assuming have navigated to the working directory in the DOS window.
% Note that all output that matlab normally sends to its command window
% now goes into the logfile run.log. If something does not work, this is
% where to look.
%
% run.m does the following
% it sets the paths for matlab, so it can fined the mfLab mfiles
% runs mf_setup.m
% exits
%
% To immediately run mf_analyze after mf_setup you must set the variable
%AFTERMFSETUP='mf_analyze'
% in mf_adapt.m, see comments in that file
%
% For intstructions how to generate the necessary files for PEST
% see mfiles
% genpestfiles.m
%
% TO 100708

% Print date
fprintf('Running mf_setup from the windows command line\n');
fprintf('Date/time = %s\n',datestr(now));

% Then sets the mfLab paths
P='Z:\tolsthoorn On My Mac\GRWMODELS\mflab\mfiles\';
addpath([P 'read']);
addpath([P 'write']);
addpath([P 'etc']);
addpath([P 'gridcoords']);
addpath([P 'fdm']);
addpath([P 'nhi2agv'])

% Then runs mf_setup.m
fprintf('Running mf_setup\n');
mf_setup

% And finally exits Matlab
fprintf('Exiting matlab, bye bye\n');
exit

% Look in the run.log file for the output of matlab that usually goes to
% the command window.