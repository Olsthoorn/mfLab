function copyFig(varargin)
% export current figure to file
% USAGE copyFig([fName][,'r',N][,'d',fmt])
% defaults fName = 'copiedFig',
%          resolution = 300 dpi
%          fmt        = 'png'
% EXAMPLES
%       copyFig();
%       copyFig(fName);     % default 
%       copyFig('-r',300);             % default fmt = 'png'
%       copyFig('-r',300,'-d','png');
%       copyFig('-r',300,'-d','TIFF',fileName);
%
% TO 141130
defName = 'copiedFig';

[res  ,varargin] = getProp(varargin,'-r',300);
[fmt  ,varargin] = getProp(varargin,{'-d','-f'},'png');
[fName,  ~ ]     = getNext(varargin,'char','');

d = dir();
id = strmatchi('Doc',{d.name});
if id && d(id).isdir,  P = d(id).name; else  P = ''; end

if isempty(fName)
    d = dir(fullfile(P,[defName '*']));
    if ~isempty(d)
        d=d(end);
        fName = d.name(1:strfind(d.name,'.')-1);
        nr    = fName(numel(defName)+1:end);
        if isempty(nr)
            nr='001';
        else
            nr = sprintf('%03d',str2double(nr)+1);
        end
    else
        nr = '001';
    end
    fName = [defName nr];
end

figNum = sprintf('-f%d',get(gcf,'Number'));
res    = sprintf('-r%d',res);
fmt    = sprintf('-d%s',fmt);
fName  = fullfile(P,fName);

print(figNum,res,fmt,fName);

fprintf('print %s %s %s %s\n',figNum,res,fmt,fName);
