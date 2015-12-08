function [value,varargin] = getThis(varargin,varClass,default)
%GETTHIS gets value of input argument of ginven class. Class must be unique in varargin.
%
% Example:
%    [value,varargin] = getThis(varargin,class,default)
%
% If not then the last is retured.
% Useful for wellObj, gridObj etc.
%
% See also: getNext getProp
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
