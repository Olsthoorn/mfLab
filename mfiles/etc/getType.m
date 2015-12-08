function [value,varargin] = getType(varargin,varClass,default)
%GETNEXT get variable whose type is varClass from varargin, if not present use default
%
% USAGE:
%    [value,varargin] = getType(varargin,class,default)
%
% the first variable in varargin whose type matches class is extracted from the
% varargin and sent to value. This variable is removed from varargin.
% used to process input arguments
%
% SEE ALSO: getNext, getProp, getType, getWord
%
% Example:
%  [ax,varargin] = getNext(varargin,'axis',gca);
%  [I ,varargin] = getnext(varargin,'gridObj',[]);
%
% varargin in output is without prop and its argument
%
% TO 130328 130807

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if isempty(varargin)
    value = default;
    return;
end

if strcmpi(varClass,'axis');   
    a = cellfun(@isaxis,varargin,'UniformOutput',false);
    for i=1:numel(a)
        if all(a{i})
            value = varargin{i};
            varargin(i)=[];
            return;
        end
    end
end

varClasses = cellfun(@class,varargin,'UniformOutput',false);

i = strmatchi(varClass,varClasses);
if ~i
    value = default;
else
    value = varargin{i(1)};
    varargin(i(1))=[];
end
