function I=strmatchi(str,strs,varargin)
%STRMATCHI finds matches of (cell array of) string(s) in cell of array of strings
%
% USAGE:
%   I=strmatchi(str  ,strss[,flag[,flag2]])  --- find location of str in strs
%   I=strmatchi(strs1,strs2[,flag[,flag2]])  --- find location of str in strs
%
% where str is a string and strs1 and strs2 are cell arrays with strings.
% 
% Almost the same as function strmatch, however base case insensitive.
% Flag:
%   'exact' to make sure the target strings match the
%           requested string exactly (but not case sensitive).
%   'once' to guaranty only one outcome (the first found)
%   'empty' to issue empty for no-hit instead of false (0)
%   'skip'  neglects hits that are NaN. If first string is a cellArray,
%           then it throws out values that are NaN.
%   A no-hit yields false instead of empty.
% if the last flag whatever then number of input argument is the
% word 'empty' then a no-hit will yield an empty matrix instead of false.
%
% example:
%    strmatchi('john',{'Johnson','and','Partners','house of John'})
%    strmatchi('john',{'Johnson','and','Partners','house of John'},'exact')
%    strmatchi('john',{'Johnson','johnny','and john''s doc','house of John'})
%
%  
%

% See also: strcmp strcmpi findstr strfind ismember
%
% TO 070630 081231 100120 110512 110807 (strs1 as input now possible)
% TO 131021 (getWord and getNext added for more structured input)

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if nargin<2, error('%s: two arguments required'); end

[exact,varargin] = getWord(varargin,'exact');
[once ,varargin] = getWord(varargin,'once');
[empty,   ~    ] = getWord(varargin,'empty');

% To allow looping in a uniform way
if ischar(str),  str = {str};  end
if ischar(strs), strs= {strs}; end

I = false(size(strs));

for i=1:numel(str)
    if exact
        for j=1:numel(strs)
            I(j) = I(j) | strcmpi(str{i},strs{j});
        end
    else
        I = I | strncmpi(str{i},strs,numel(str{i}));
    end
    if once && any(I)
        break;
    end
end

I=find(I);

if isempty(I)
    if ~empty
        I = false;
    end
    return;
else
    if once
        I=I(1);
    end
end

