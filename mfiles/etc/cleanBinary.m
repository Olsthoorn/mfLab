function cleanBinary(fname,varargin)
%CLEANBINARY  removes record length info from binary files
%
% USAGE:
%    cleanBinary(fname[,verbose])
%
% Removes record lenght bytes from binary if present
%    This is requried if a binary file was generated using
%    ACCESS/'SEQUENTIAL'/ instead of /'STREAM'/ while the reading file
%    assumes /'STREAM'/ as input.
%
%    it used in mf_setup to clean the binary heads and budget file before
%    MPATH is run, as otherwise it would crash if record length information
%    is in these files. Whether this information is put in the file depends
%    on the fortran compiler used to compile the executables. The fortran
%    compiler used lies outside the power of mflab, it is given, so
%    cleanBinary is a way to handle it.
%
% TO 120218

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

verbose = nargin>1;

if verbose
    fprintf('Cleanding binary file %s, if necessary ...\n',fname);
end

fpi = fopen(fname,'r');

fseek(fpi,0,1); LFile = ftell(fpi); frewind(fpi); % end of file

n=0;
while ftell(fpi)<LFile && n<5
    I1 = fread(fpi,1,'int32');
    n=n+1;
    fread(fpi,I1,'uchar');
    I2 = fread(fpi,1,'int32');
    if I1~=I2,
        if verbose
            fprintf('File <<%s>> has no reclen info, it is streamed file.\n',fname);
            fseek(fpi,0,1);
            fprintf('Number of bytes <<%d>>.\n',ftell(fpi));
        end
        fclose(fpi);
        return;
    end
end

fseek(fpi,0,1);

%% File info
if verbose
    fprintf('File <<%s>> contains record length bytes, which will be removed ...\n',fname);
    fprintf('File length in bytes <<%d>>.\n',ftell(fpi));
end

fclose(fpi);

%% Saving old file
fnameOld = [fname '.old'];
if verbose, fprintf('Renaming file <<%s>> to <<%s>>\n',fname,fnameOld); end
movefile(fname,fnameOld);

%% Opening files
fpi=fopen(fnameOld,'r');
fpo=fopen(fname   ,'w');

%% Processing
n=0;
while ftell(fpi)<LFile
    I1 = fread(fpi,1,'int32');
    n=n+1;
    fwrite(fpo,fread(fpi,I1,'uchar'),'uchar');
    I2 = fread(fpi,1,'int32');
    if I1~=I2,
        error('%s: binary file %s, record %d, ftell=%d (length at beginning %d ~= length at and %d\n',...
            mfilename,fname,n,ftell(fpi),I1,I2);
    end
end

%% File info
if verbose
    fprintf('file <<%s>> has its record length bytes removed.\n',fname);
    fprintf('Number of records in file <<%d>>.\n',n);
    fprintf('File length in bytes <<%d>>.\n',ftell(fpo));
end

%% Cloaing files
fclose(fpo);
fclose(fpi);

% finally remove the orginal file
delete(fnameOld);


