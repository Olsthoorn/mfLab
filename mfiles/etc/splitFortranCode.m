function [n,code,w,d] = splitFortranFormat(ftfmt)
%SPLITFORTRANCODE -- split up the fortran code in [n,code,width,decimals]
%
%USAGE:
%   [n,code,w,d] = splitFortranCode(ftfmt)
%
% examples
%    [n,code,w,d] = splitFortranCode('(20F10.3)');
%    [n,code,w,d] = splitFortranCode('(20F10.3)');
%    [n,code,w,d] = splitFortranCode('(20F.3)');
%    [n,code,w,d] = splitFortranCode('(20e10.3)');
%    [n,code,w,d] = splotFortranCode)'(10I2)');
%
% remove blanks and brackets
ftfmt = ftfmt(ftfmt~=' ' & ftftm ~= '(' & ftfmt ~= ')');

iCode = regexpi('[idefg]',ftfmt);
iDot  = regexpi('\.',ftfmt);
n     = sscanf(ftfmt(1:iCode-1),'%d',1);
code  = ftfmt(iCode);
w     = sprintf('%d',ftfmt(iCode+1,:));
d     = sprintf('%d',ftfmt(iDot +1,:));
