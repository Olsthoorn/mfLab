function [value,varargin] = getNext(varargin,varClass,default)
%GETNEXT get prop value from varargin, if not present use default
%
% USAGE:
%    [value,varargin] = getNext(varargin,class,default)
%
% the first argument of varargin, varargin{1} is put into value if it is of
% the requrested class. If true, then then varargin(1) is removed. If
% false varargin is left untouched and value becomes the given default.
% Used often to process input arguments
%
% SEE ALSO: getNext, getProp, getType, getWord
%
% Example:
%  [ax,varargin] = getNext(varargin,'axis',gca);
%  [I ,varargin] = getnext(varargin,'double',[]);
%
% varargin in output is without prop and its argument
%
% TO 130328

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if isempty(varargin)
    value = default;
elseif strncmp(varClass,'axis',numel(varClass)) % special, case not a class
    if isaxis(varargin{1})
        value = varargin{1};
        varargin(1)=[];
    else
        value = default;
    end    
elseif ~any(strmatchi(class(varargin{1}),varClass))
    value = default;
else
    value = varargin{1};
    varargin(1)=[];
end
