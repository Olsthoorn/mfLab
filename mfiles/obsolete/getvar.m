function [value,varargin] = getvar(prop,varargin,default)
% [value,varargin] = getvar(prop,varargin,default)
% get prop value from varargin, if not present use default
% varargin in output is without prop and its argument
%
% TO 121219

i = strmatchi(prop,varargin,'exact');
if i
    value = varargin{i+1};
    varargin(i:i+1)=[];
else
    value = default;
end
