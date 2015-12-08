%UNPACK script to to get variables from ModelObj array into the workspace (through eval)
%
% Example:
%    unpack
%
% The model is an array of modelObj with fields name, var, type
% Each element holds one object or array
% see modelObj
%
% See also: modelObj
%
% TO 120813

for i=1:numel(Model)
    fprintf('Putting    %-20s  of mfLab type     %-15s   in the workspace\n',Model(i).name,Model(i).type);
    if ~isempty(Model(i).name)
        eval([ Model(i).name '= Model(i).var;']);
    end
end
