function unpack(o,varargin)
%UNPACK script to to get variables from ModelObj array into the workspace (through eval)
%
% USAGE:
%    Model.unpack
%    Model.unpack('type','gridObj');
%    Model.unpack('name','gr')
%    Model.unpack('type',{'gridObj','well'},'name',{'HK','VK'})
%
% The model is an array of modelObj with fields name, var, type
% Each element holds one object or array
% see modelObj
%
% See also: modelObj
%
% TO 120813 131120

[type,varargin] = getProp(varargin,'type','');
[name,varargin] = getProp(varargin,'name','');

if ~isempty(varargin)
    fprintf('varargin not completely used');
end

if ~isempty(type)
    o = o(ismember({o.type},type));
end
if ~isempty(name)
    o =o(ismember({o.name},name));
end

for io=1:numel(o)
    fprintf('Putting    %-20s  of mfLab type     %-15s   in the workspace\n',o(io).name,o(io).type);
    if ~isempty(o(io).name)
        eval([ o(io).name '=  o(io).var ;' ]);
    end
end

save2base(true);
