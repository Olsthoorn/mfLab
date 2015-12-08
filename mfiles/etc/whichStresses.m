function strss = whichStresses(varargin)
%WHICHSTRESSES -- finds which stress packages have been used in current simulation
%
% USAGE:
%   strss = whichStresses()
%   strss = whichStersses(stresses)
%
%   stresses is a cell array of stresses to look for if omitted the default
%   list is used as defined below.
%
%   strss is a cell array with the names of the stress packages that have
%   been used in the current simulation. The possible stress packages are
%   RCH EVT ETS WEL MNW1 MNW2 DRN RIV GHB HFP
%
%   is find the most recent .man file and looks for the stress packages
%   defined inside and decides whether they are one of the given lies of
%   stress packages
%
%   TO 151126

if nargin<1
    legalStresses = {'CHD' 'RCH' 'EVT' 'ETS' 'WEL' 'MNW1' 'MNW2' 'DRN' 'RIV' 'GHB'};
else
    legalStresses = varargin{1};
end

%% Get the most recent nam file (this was run)

d = dir('*.nam'); % get all .nam files

if isempty(d)
    error('No .nam files in this directory, can''t continue');
end

% only keep the most recent one
[~,I] = sort([d.datenum]);
d     = d(I(end));

% scan the file
fid = fopen(d.name,'r');

strss = false(size(legalStresses));
while 1
    l      = fgets(fid);  if l==-1, break; end

    pkg    = l(1:(strfind(l,' ')-1));
    i      =  find(ismember(legalStresses,pkg));
    if ~isempty(i)
        strss(i) = true;
    end
end

strss = legalStresses(strss);
