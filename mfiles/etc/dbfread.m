function data=dbfread(dbffilename,varargin)
%DBFREAD reads data from dbf file
%
% USAGE:
%    data=dbfread(dbffilename)
%    data=dbfread(dbffilename [,'v'])
%
%    dbffilename is the completepath\filename needed with or without '.dbf'.
%    the oprion 'v' means verbose, i.e. more info about contents of the dbf
%    file. Default is non-verbose.
%
%    data(NField) is a structarray holding the contents of the dbffile
%      data(iField) contains the fieldnames and meta values
%      data(iField).values contains the record values of the field
%      data(iField).values is a numeric array if field type is numeric
%          otherwise it is a cell array.
%
%    The dbf specification is given at the end of this m-file
%
%  TO 970901 FS 000101 TO 090729 - complete overhaul

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later
verbose = getWord(varargin,{'v','verbose'});

if verbose, fprintf('\n\ndbfread --- %s\n',datestr(now)); end

if length(dbffilename)>4 && strcmp(dbffilename(end-3:end),'.dbf')    % don't worry about '.dbf'
    dbffilename=dbffilename(1:end-4);
end

fp =fopen([dbffilename,'.dbf'],'r');
if fp<0, error('ERROR in DBFREAD: Can''t open file ''%s.dbf''\7\n',dbffilename); end

% reading dbf file bytewise according to the meaning of the bytes
           fread(fp,1,'int8');   % 1 byte valid dBASE III Plus table file 03h without and 83h with a memo .dbt file
yymmdd   = fread(fp,3,'int8');   % [yy mm dd]' file date
NRec     = fread(fp,1,'int32');  % number of records in dbf file
begindat = fread(fp,1,'int16');  % start position of first data record (in bytes, i.e. first=0 not 1)
reclen   = fread(fp,1,'int16');  % data record length

if verbose
    fprintf('Reading dbf file ''%s'',\nwhich was last updated on %s\n',dbffilename,...
        datestr(datenum(1900+yymmdd(1),yymmdd(2),yymmdd(3),0,0,0)));
end
%% Read subreads (field data)

if verbose
    fprintf('%10s  %-5s  %s\n','fieldlen','fieldtype','fieldname');
end

NField=(begindat-1-32)/32;
data(NField).fieldname='dummy to allocate data array';
for iField=1:NField
    fseek(fp,iField*32,'bof');
    data(iField).fieldname           = strtrim(fread(fp,[1,11],'*char'));  % first 11 bytes of subrecord
    data(iField).fieldtype           = fread(fp,1,'*char'); % field type, 1 byte
    data(iField).fielddisplacement   = fread(fp,1,'int32'); % displacement of field inrecord
    data(iField).fieldlen            = fread(fp,1,'uint8'); % field length in bytes
    data(iField).decplaces           = fread(fp,1,'uint8'); % number of decimal places
    data(iField).fieldflags          = fread(fp,1,'uint8'); % field flags
    data(iField).autoincrementvalue  = fread(fp,1,'int32'); % autoincrement next value
    data(iField).autoincrementstepval= fread(fp,1,'int8' ); % value of autoincrement step
    switch data(iField).fieldtype
        case {'D', 'T'}
            data(iField).values=NaN(NRec,3);
        case {'Y', 'N', 'F', 'B', 'I', 'L'}
            data(iField).values=NaN(NRec,1);
        case {'C', 'M', 'G', 'P'}
            data(iField).values=cell(NRec,1);
        otherwise
            error('Unknown field type ''%s'', in field %s, fieldnr %d\n',...
                data(iField).fieldtype,data(iField).fieldname,iField);
    end
    if verbose
        fprintf('%10d  %-5s  %s\n',data(iField).fieldlen,data(iField).fieldtype,data(iField).fieldname)
    end
end


% count reclen form fieldlen of individual fields
% FL=0;
% for iField=1:NField
%     FL=FL+data(iField).fieldlen;
%     fprintf(' %d',data(iField).fieldlen)
% end
% fprintf('\nTotal record length %d\n',FL);

%% Read the records
for iRec=1:NRec
    startrecdata=begindat+(iRec-1)*reclen+1;      % compute file pointer position of this record
    fseek(fp,startrecdata,'bof');                 % position file pointer
    for iField=1:NField                           % run over the fiels using the fieldlen
        s=char(fread(fp,[1,data(iField).fieldlen],'uchar')');  % store field bytes for interpretaion
        if any(deblank(s))  % s must not be empty
            switch data(iField).fieldtype
                case 'C'   %   Character
                   data(iField).values(iRec)={s'};
                case 'Y'   %   Currency
                   data(iField).values(iRec)=sscanf(s,'%f',1);
                case 'N'   %   Numeric
                   data(iField).values(iRec)=sscanf(s,'%f',1);
                case 'F'   %   Float
                   data(iField).values(iRec)=sscanf(s,'%f',1);
                case 'D'   %   Date
                   data(iField).values(iRec,:)=sscanf(s,'%d',[1,3]);
                case 'T'   %   DateTime
                   data(iField).values(iRec,:)=sscanf(s,'%s',[1,3]);
                case 'B'   %   Double
                   data(iField).values(iRec)=sscanf(s,'%f',1);
                case 'I'   %   Integer
                   data(iField).values(iRec)=sscanf(s,'%d',1);
                case 'L'   %   Logical
                   data(iField).values(iRec)={s};
                case 'M'   %   Memo
                   data(iField).values(iRec)={s};
                case 'G'   %   General (bytes, use %char, not %*char)
                   data(iField).values(iRec)={s};
                case 'P'   %   Picture  (bytes, use %char, not %*char)
                   data(iField).values(iRec)={s};
                otherwise
                    error('Unknown field type ''%s'' of field ''%s'' %d',data(iField).fieldtype,data(iField).fieldname,iField);
            end
        end
    end
end
        
fclose(fp);

% DBF File structure from http://www.dbf2002.com/dbf-file-format.html (date 090729)
% A DBF file consists of a header record and data records. The header record defines the structure of the table and contains any other information related to the table. The header record starts at file position zero. Data records follow the header, in consecutive bytes, and contain the actual text of the fields. 
% 
% Note   The data in the data file starts at the position indicated in bytes 8 to 9 of the header record. Data records begin with a delete flag byte. If this byte is an ASCII space (0x20), the record is not deleted. If the first byte is an asterisk (0x2A), the record is deleted. The data from the fields named in the field subrecords follows the delete flag.
% The length of a record, in bytes, is determined by summing the defined lengths of all fields. Integers in table files are stored with the least significant byte first.
% 
% DBF File Header
% Byte offset Description 
% 0 File type:
% 0x02   FoxBASE
% 0x03   FoxBASE+/Dbase III plus, no memo
% 0x30   Visual FoxPro
% 0x31   Visual FoxPro, autoincrement enabled
% 0x32   Visual FoxPro with field type Varchar or Varbinary
% 0x43   dBASE IV SQL table files, no memo
% 0x63   dBASE IV SQL system files, no memo
% 0x83   FoxBASE+/dBASE III PLUS, with memo
% 0x8B   dBASE IV with memo
% 0xCB   dBASE IV SQL table files, with memo
% 0xF5   FoxPro 2.x (or earlier) with memo
% 0xE5   HiPer-Six format with SMT memo file
% 0xFB   FoxBASE 
% 1 - 3 Last update (YYMMDD) 
% 4  7 Number of records in file 
% 8  9 Position of first data record 
% 10  11 Length of one data record, including delete flag 
% 12  27 Reserved 
% 28 Table flags:
% 0x01   file has a structural .cdx
% 0x02   file has a Memo field
% 0x04   file is a database (.dbc)
% This byte can contain the sum of any of the above values. For example, the value 0x03 indicates the table has a structural .cdx and a Memo field. 
% 29 Code page mark 
% 30  31 Reserved, contains 0x00 
% 32  n Field subrecords
% The number of fields determines the number of field subrecords. One field subrecord exists for each field in the table. 
% n+1 Header record terminator (0x0D) 
% n+2 to n+264 A 263-byte range that contains the backlink, which is the relative path of an associated database (.dbc) file, information. If the first byte is 0x00, the file is not associated with a database. Therefore, database files always contain 0x00. 
% 
% Field Subrecords Structure
% Byte offset Description 
% 0  10  Field name with a maximum of 10 characters. If less than 10, it is padded with null characters (0x00). 
% 11 Field type: 
% C      Character
% Y      Currency
% N      Numeric
% F      Float
% D      Date
% T      DateTime
% B      Double
% I      Integer
% L      Logical
% M      Memo
% G      General
% C      Character (binary)
% M      Memo (binary)
% P      Picture 
% 12  15 Displacement of field in record 
% 16 Length of field (in bytes) 
% 17 Number of decimal places 
% 18 Field flags:
% 0x01   System Column (not visible to user)
% 0x02   Column can store null values
% 0x04   Binary column (for CHAR and MEMO only) 
% 0x06   (0x02+0x04) When a field is NULL and binary (Integer, Currency, and Character/Memo fields)
% 0x0C   Column is autoincrementing 
% 19 - 22 Value of autoincrement Next value 
% 23 Value of autoincrement Step value 
% 24  31 Reserved 
% 
