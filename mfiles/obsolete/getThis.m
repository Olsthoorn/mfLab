function [value,varargin] = getThis(varargin,varClass,default)
% [value,varargin] = getThis(varargin,class,default)
% get value of for class given. Class must be unique in varargin.
% If not then the last is retured.
% Useful for wellObj, gridObj etc.
%
% TO 130328

for i=numel(varargin):-1:1
    if strcmpi(class(varargin{i}),varClass)
        value = varargin{i};
        varargin(i)=[];
        return;
    end
end

value = default;
