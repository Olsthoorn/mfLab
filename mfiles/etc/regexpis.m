function L = regexpis(cellArrayOfStr,exp,varargin)
%REGEXPIS find strings in cellarrayOfStr that full file the regular expression exp.
%
% USAGE
%    L = regexpis(cellArrayOfStr,exp)
%
% Example
%    L = regexpis({'john','mary',george'},j*|g*,once);
%
% See also: regexp regexpi
%
% TO 120507

if ~iscell(cellArrayOfStr)
    error('%s: first arg must be a cell array of strings',mfilename);
end
if ~ischar(exp)
    error('%s: second arg must be a regular expression',mfilename);
end

L = find(~cellfun(@isempty,regexpi(cellArrayOfStr,exp,'once'))); %,'once'),varargin{:});
