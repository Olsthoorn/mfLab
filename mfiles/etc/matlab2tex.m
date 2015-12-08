function status = matlab2tex(mfName,varargin)
%MATLAB2TEX -- write TEX files for matlab m-file
% USAGE:
% function status = matlab2tex(mfile,'lineNumbers')
%  Writes a file with tex for syntax colored listing of matlab code in mfile
%  an example of how to include tex in Lyx: 
%  Paste into Document preamble (without matlab comments, of course)
%      \usepackage{alltt}
%      \usepackage{color}
%      \definecolor{string}{rgb}{0.7,0.0,0.0}
%      \definecolor{comment}{rgb}{0.13,0.54,0.13}
%      \definecolor{keyword}{rgb}{0.0,0.0,1.0}
%  Insert a float figure
%  Click in figure, start an insert ERT (ctrl-L)
%  Do Insert|file|plain text and navigate to tex file
%  click outside ERT
%  ---
%  inputs: 
%    mfile: the m-file path and mfName. If no path is specified, looks in matlab working dir 
%    linenumbers: (optional = false) If present, includes linenumbers in listing
%  outputs:
%    creates a file mfilemfName.tex in same dir as original matlabfile
%    status: the last 'status' from the file open/close commands. 0 if executed without errors
%  REA 5/11/09

% options for highlight function
 opt.type = 'tex';
[opt.linenb, varargin] = getWord(varargin,'line');


[path, mfName, ext] = fileparts(mfName); % ext ignored, '.m' assumed

if ~isempty(ext)
    if ~strcmp(ext,'.m'); error('no m file given'); end
end

texFname =  [path 'tex' filesep mfName '.tex'];

mfName = [path mfName '.m'];

if isempty(varargin)
    if ~exist([path 'tex'],'dir')
        mkdir([path 'tex']);
    end

    fidOut = fopen(texFname,'wt');
    if fidOut>0
        highlight(mfName,opt,fidOut);
    else
        error('Can''t generate tex file\n');
    end
    status = fclose(fidOut);
    fprintf('%s generated\n',texFname);
else
    if ~isnumeric(varargin{1})
        error('last argument should be numeric, which implies printing to screen');
    end
    highlight(mfName,opt,1);
end

function [b,varargin] = getWord(varargin,word)
    i = strncmpi(word,varargin,2);
    b = i>0;
    if b
        varargin(i)=[];
    end

