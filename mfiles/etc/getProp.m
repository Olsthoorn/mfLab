function [value,varargin] = getProp(varargin,propertyName,default)
%GETPROP get property value from varargin, if not present use default
%
% USAGE:
%    [value,varargin] = getProp(varargin,propertyName,default)
%
% often used to make processing input robustly. This function
% grabs the propertyName,value from varargin. It grabs the
% propertyName,value pair from anywhere in the varargin and then removes
% them from varargin.
%
% if string propertyName is found in varargin{:}, the next value, i.e.
% varargin{i+1} is transferred to the output and varargin(i:i+1) are
% removed. If the propertyName does not exist in varargin, then varargin is
% untouched, while value becomes default.
%
% SEE ALSO: getNext, getProp, getType, getWord

% varargin in output is without prop and its argument
%
% TO 130328

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if isempty(varargin)
    value = default;
else
    i = strmatchi(propertyName,varargin);
    if numel(i)>1
        i = strmatchi(propertyName,varargin,'exact');
    end
    if numel(i)>1
        error('%s: Property names not unique %s',mfilename,sprintfs(' <<%s>>',varargin{i}));
    end
    if i
        value = varargin{i+1};
        varargin(i:i+1) = [];
    else
        value = default;
    end
end
