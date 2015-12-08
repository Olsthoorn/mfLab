function A = splitFortranFormat(ftfmt)
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
% where code is the ascii value of the leter denoting the format
% code, i.e. the fortran format code is char(A(2))
%
% 140113

% remove blanks and brackets
ftfmt = ftfmt(ftfmt~=' ' & ftfmt ~= '(' & ftfmt ~= ')');

A = sscanf(ftfmt,'%d %c %d %c 5d');
if A>4, A(4)=A(5); end

