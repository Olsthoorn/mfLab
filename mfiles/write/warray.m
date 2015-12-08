function warray(fid,var,unit,ftfmt,txt,ctrlRecDesired,freeFmt,FName)
%WARRAY writes MODFLOW and MT3DMS arrays to modflow input files (general array writer)
% General array writer. Does what UD2REL UD2INT UD1REL and UD1INT does in Modflow.
%
% Example:
%    warray(fid,var,unit,ftfmt,txt[,ctrlRecDesired?[,freeFmt?[,FName]]]) --- write out a 1D or 2D MT3D matrix
%
% ctrlRecDesired?  if omitted, empty or true requires control record (default)
% freeFmt?     if omitted, empty or false uses fixed format, if true uses freeFmt format
% FName     required if data has to be read from separate file using OPEN/CLOSE
%
% Example:
%  warray(fid,var,unit,ftfmt,txt,true,false);  % --- fixed format records and header txt is optional info string
%  warray(fid,var,unit,ftfmt,txt ,true  );  % --- freeFmt format using  EXTERNAL
%   instead of INTERNAL we use EXTERNAL with an already open unit
%  warray(fid,var,unit,ftfmt,ext,true,Fname); % --- freeFmt format using OPEN/CLOSE
%
% The freeFmt format option always used F15.7 format for floating point numbers, max resolution
%
% ftfmt is a fortran format (10E3.2) (4I1) etc  10G6.1 also works
%
% TO 070703 081228 090907 100911 130206 132928
%
% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under freeFmt software foundation GNU license version 3 or later
%
% control record fortran format string, replace g by e or
% otherwise MODFLOW won't understand it. Leave ftfmt intact
% because it must be used as is in foldprint. But wel will reduce the
% number of values to print on a line to the mininum of size(var,2) and the
% number specified in the Fortran format.

    if nargin<6, ctrlRecDesired=true;  end
    if nargin<7, freeFmt=false;        end;     if isempty(freeFmt), freeFmt = false; end
    
    if ~islogical(ctrlRecDesired), error('%s: ctrlRecDesired must be of class logical',mfilename); end
    if ~islogical(freeFmt),        error('%s: freeFmt must be of class logical',mfilename); end
    

    IPRN = 3;

%    ftfmtstr=ftfmt; ftfmtstr(lower(ftfmtstr)=='g')='e';

    [~,~,ftfmt]=ftfmt2ml(ftfmt,size(var,2));

    if nargin<5 || isempty(txt), txt=''; end

    var(isnan(var))=0; % USGS models will not accept NaN in input

    if ~ctrlRecDesired
        % just print the data and return
        foldprint(fid,var,ftfmt,txt);
        return;
    end

    %% Control record format and string

    % Are all variables constant ??
    constant = all(var(:)==var(1));
    if constant
        unit  = 0;
        CONST = var(1);
    else
        CONST = 0;
    end

    if ~freeFmt
        if ~constant
            % old fashioned control record
            fprintf(fid,'%10d%10d  %-18s%10d     %s\n',unit,CONST,ftfmt,IPRN,txt);
            % followed by the data
            foldprint(fid,var,ftfmt,txt);
        else
            IPRN=-1;
            % does CONST fit within 10 characters?
            s=sprintf('%g',CONST);  % see what would be printed
            if length(s)>10         % if more than 10 characters
                i=strfind(s,'e');   % find end of field before exponent
                if ~isempty(i)
                    % truncate and assemble to field of 10 char long
                    s=[s(1:(10-(length(s)-i+1))) s(i:end)];
                else
                    s=s(1:9);
                end
            end
            % Just print the control record:
            fprintf(fid,'%10d%10s  %-18s%10d     %s\n',...
                unit,s,ftfmt,IPRN,txt);   % print with float as string of 10 char wide
        end
    else  % free format printing
        if constant
            fprintf(fid,'CONSTANT %12g     %s\n',CONST,txt);
        else
            if nargin<8
                if unit==0  % use INTERNAL
                    % print control record without unit number
                    fprintf(fid,'INTERNAL %12g %s %d     %s\n',CONST,ftfmt,IPRN,txt);
                else        % use EXTERNAL
                    fprintf(fid,'EXTERNAL %d %12g %s %d     %s\n',unit,CONST,ftfmt,IPRN,txt);
                end
                % print the data
                foldprint(fid,var,ftfmt,txt);
            else  %use 'OPEN/CLOSE'
                % print control record
                fprintf(fid,'OPEN/CLOSE %s %12g %s %d     %s\n',FName,CONST,ftfmt,IPRN,txt);
                % print the data
                foldprint(fid,var,ftfmt,txt);
            end
        end
    end
end

function foldprint(fid,var,ftfmt,txt)
    % foldprint(fid,var,ftfmt [,norec [, txt])
    % prints according to fortran format ftfmt
    % if norec is given (arbitrary) var is printed as a list else it is printed
    % as an array
    % if text is present, the text is printed at the end of the list.
    % row by row with a newline behind every row.

    [mlfmt NW]=ftfmt2ml(ftfmt,[]);
    
    if mlfmt(end)=='i', mlfmt(end)='d'; end
    
    [Nrow,Ncol]=size(var);
    
    NW = min(NW,Ncol);

    for iRow= 1:Nrow
        k1=1;
        while k1<=Ncol
            k2=min(k1+NW-1,Ncol);
            fprintf(fid,mlfmt,var(iRow,k1:k2));
            if k2==Ncol
%                fprintf(fid,'    %s\n',txt);
                fprintf(fid,'\n');
            else
                fprintf(fid,'\n');
            end
            k1=k1+NW;
            if k2>=Ncol
                break;
            end
        end
    end
end
