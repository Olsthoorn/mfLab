function save2base(overwrite)
%
% save_to_base() copies all variables in the calling function to the base
% workspace.  This makes it possible to examine a function's internal
% variables from the Matlab command prompt after the calling function
% terminates.
%
% If the optional argument "overwrite" is present and has a non-zero value,
% variables in the workspace of the caller will overwrite variables in the
% base workspace having the same name.  Otherwise, preexisting variables in
% the base workspace will not be overwritten.
%
% Phillip M. Feldman

% Version 1.1, 18 Aug. 2009: Added fix to exclude 'ans' variable.
%
% Version 1.0, 2 June 2009: Initial version.
%
% Note: I'd like to acknowledge assistance from Sean Little of the
% MathWorks.

% get variable names in context of caller of current function
ws_caller= evalin('caller','who()');

if nargin>0 && overwrite
   % Overwrite preexisting variables:
   variables= ws_caller;
else
   % Do not overwrite preexisting variables:
   % get variable names in basic (Matlab's) workspace
   ws_base= evalin('base','who()');
   
   % get names that are in the caller but not in the workspace
   variables= setdiff(ws_caller, ws_base);
end

for i= 1 : length(variables)
   if ~strcmp(variables{i},'ans') % skip 'ans'
      tmp= evalin('caller',variables{i});
      
      % assign variables{i} in caller to name tmp in workspace
     
      assignin('base'  ,variables{i},tmp);
   end
end
