function ftfmt=mlfmt2ft(mlfmt,N)
%MLFMT2FT converts a matlab number format to a fortran equivalent.
%
% Example:
%    ftfmt=mlfmt2ft(mlfmt,N);
%    ftfmt=mlfmt2ft('%12.4f',8);
%    yields '(8F12.4)'
%
% TO 081225

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if isempty(mlfmt), error('format not specified!'); end

if isempty(findstr('%',mlfmt)),
    error(['Illegal matlabformat: ',mlfmt]);
else
    mlfmt=mlfmt(findstr('%',mlfmt):end);  % remove preceeding blanks
end

A=sscanf(mlfmt,'%c %d %c %d %c %d %c');
if A(1)~='%', error(['Illegal natlab format: ',mlfmt]); end

switch length(A)
    case {1, 2}
        if mlfmt(end)=='d'
            ftfmt='(10I10)';
        else
            ftfmt='(10E15.5)';
        end
    case 3
        ftfmt=sprintf('(%d%c%d)'    ,N,upper(char(A(3))),A(2));
    case 4
        ftfmt=sprintf('(%d%c%d)'    ,N,upper(char(A(4))),A(2),A(3));
    otherwise
        if A(3)~='.', message(['Illegal matlab format: ',mlfmt]); end
        ftfmt=sprintf('(%d%c%d%c%d)',N,upper(char(A(5))),A(2),A(3),A(4));
end
