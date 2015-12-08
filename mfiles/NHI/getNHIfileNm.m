function [fname, names]= getNHIfileNm(xlsname,sheetnm,columnHdr,VarNm)
%GETNHIFILENM gets NHI file name for workbook saved in mfLab
%
% Example:
%    fname = getNHIfileNm((xlsname,sheetnam,filesdir,VarNm) -- retrieves the file name for this NHI variable
%
% This function is used for variables that are in only one file such as
% WELL or for any file with an exact variable size of NHI or of a submodel
% the variable of which is entirely contained in this file.
%
% TO 120427

warning('off') %#ok
[~,~,txt] = xlsread(xlsname,sheetnm,'','basic');
warning('on'); %#ok
hdr = txt(1,    :);
txt = txt(2:end,:);

jVar         = strmatchi(columnHdr,hdr);
jZip         = strmatchi('zip'  ,hdr);
jFile        = strmatchi('file' ,hdr);
jDescription = strmatchi('descr',hdr);

I= strmatchi(VarNm,txt(:,jVar),'exact');   % col jVar contains file names

names={txt{I,jVar},txt{I,jFile},txt{I,jZip},txt{I,jDescription}};

if I==0,
    error('getNHIfileNm:scan:NotFound',....
        '%s: Can''t find file for variable<<%s>>.',mfilename,VarNm);
end

if length(I)>1
    error('getNHIfileNm:scan:FileNotUnique',...
        '%s: Files for variable <<%s>> are not unique.',mfilename,VarNm);
end

fname = txt{I,jVar+1};
    