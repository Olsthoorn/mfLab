function [mlfmt,N,ftfmt]=ftfmt2ml(ftfmt,NCOL)
%FTFMT2ML converts fortrans format to matlab format for printing
%
% Example:
%    [mlfmt,N,ftfmt]=ftfmt2ml(ftfmt[,FREE]);
%
% if FREE is true, then free format will be used for floating point
% numbers, which is e15.7 or %15.7e in Matlab
%
% Converts fortran format to matlab format
% Example:
%    [mlfmt,N,ftfmt]=ftfmt2ml('(8F7.2)');
%    yields ['%7.2f',8,ftfmt]
%
% ftfmt is adepted when neceessary to make room for the required number of
% significant digits
%
% Other ftfmt examples are (20I4) (8E15.3) (9G10.3)
%
% See also: mlfmt2ft
%
% TO 081231 120206


% Copyright -2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

ftfmt=strtrim(ftfmt); if ftfmt(1)=='(', ftfmt=ftfmt(2:end-1);  end
A=sscanf(ftfmt,'%d %c %d %c %d');
if nargin>1 && ~isempty(NCOL)
      A(1)=min(A(1),NCOL);
end

switch length(A)
    case {1 2}
        error(['Illegal fortran format: ',ftfmt]);
    case 3
%         if ismember(char(A(2)),'idID')
%             A(1) = min(A(1),40);
%         end
        ftfmt=sprintf('%d%c%d',A);
        mlfmt=sprintf('%%%d%c', A(3), lower(char(A(2))));
    otherwise  % i.e. 5 args like 10E12.5
        A(2)=lower(char(A(2)));
        if A(2)=='e' || A(2)=='g',
            A(5)=min(A(3)-6,A(5)); % TO: 6 was 7 until 120524
        end
        ftfmt=sprintf('%d%c%d%c%d',A);
        mlfmt=sprintf('%%%d.%d%c', A(3),A(5),A(2));
end

N=A(1);

ftfmt=['(' ftfmt ')'];
