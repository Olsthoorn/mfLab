function C=readMT3D2(fname,userPeriods,userLayers,userRows,userCols)
%READMT3D2 reads unformatted MT3DMS concentration file (MT3D00n.UCN) into struct 
%
% USAGE:
%   C=readMT3D([-]fname [,userPeriods [,userLayers [,userRows [,userCols]]]]]])
%   C=readMT3D(fname,-1)
%      to get info on the contents of the file.
%
% make subselections using the optional arguments
% they may be omittted as shown or used empty '' or [], for instance
%
%   C=readMT3D(fname,'',[1,3,6],'',3)
%
% fetches all transport time steps , layers 1,3 and 6, all userRows  and column 3
%
%   C=readMT3D(fname,[1 3 5:8],'',3:10,12:21}
%
% fetches transport steps 1,3 and 5..8, all layers, userRows 3..10 and
% columns 12..21
%
% put '+' before filename to get more verbose output while processing
%
% C is a cell array whose length is the number of records in the file
% it has the following fields:
%    .values is a 3D array containing the requested concentration values
%    .label        label telling the contents
%    .trpstp       transport step number
%    .period       stress period number
%    .tstp         time step number within stress period
%    .time         time since begin of simulation
%    .NCOL         number of original columns in file fname
%    .NROW         number of original Rows    in file fname
%    .NLAY         number of original layers  in file fname
%    .cols         list of fetched columns
%    .rows         list of fetched userRows
%    .lays         list of fetched layers
%
%   SEE ALSO: readbud readDat
%
% TO 090104 091214 100819
%
% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if     fname(1)=='+' || fname(1)=='1', fname=fname(2:end); verbose=1;
elseif fname(1)=='-' || fname(1)=='0', fname=fname(2:end); verbose=0;
else   verbose=0;
end

fp=fopen(fname); if fp<1, error('READMT3D: Can''t find or open file %s!',fname); end

fprintf('Reading MT3DMS binary output file <<%s>> verbose= %d\n',fname,verbose);

%% Advance a single cell-by-cell record & compute the number or records in file

bytes=contrRec(fp);  % bytes is compiler specific UNFORMATTED file record heading

nbyt=ftell(fp); fseek(fp,0,'eof'); Nbyt=ftell(fp); nRec=Nbyt/nbyt; % number of records in file

if rem(Nbyt,nbyt),
    error('UNFORMATTED file %s non standard. Sorry, can''t interpret it. Try using a better compiler.',fname);
end

%% Get the information about the contents of the file from record headings

trpStpsInFile = zeros(nRec,1);
periodsInFile = zeros(nRec,1);
tstpsInFile   = zeros(nRec,1);
labelsInFile  = cell( nRec,1);
timeInFile    = zeros(nRec,1);
layersInFile  = zeros(nRec,1);
rowsInFile    = zeros(nRec,1);
colsInFile    = zeros(nRec,1);

fprintf('Scanning headers...\n');
for i=1:nRec
    fseek(fp,(i-1)*nbyt,'bof');  % move to start of last layer in file
    
    [trpStpsInFile(i),...
        tstpsInFile(i),...
        periodsInFile(i),...
        timeInFile(i),...
        labelsInFile{i},...
        layersInFile(i),...
        rowsInFile(i),...
        colsInFile(i)...
        ]=contrRec(fp,bytes);

    if verbose
        if rem(i, 100)==0, fprintf('.'); end
        if rem(i,5000)==0, fprintf('%d records read\n',i); end
    end
end
fprintf('...finished, %d records read\n',i);

%% Evaluate what we've got in terms of the periods, timeInFile steps unique lables etcetera before we continue to select the actual data

NTRP = length(unique(trpStpsInFile));  % does not work in MT3DMS
NPER = length(unique(periodsInFile));
NLAY = length(unique(layersInFile));
NROW = max(rowsInFile);
NCOL = max(colsInFile);
TIME = unique(timeInFile);

fprintf('File contains the following:\n');
fprintf('Number of records in file : %10d\n', nRec);
fprintf('Number of transport steps : %10d\n', NTRP);
fprintf('Number of stress periods  : %10d\n', NPER);
fprintf('Highest time in file      : %10g\n', max(TIME));
fprintf('Number of times in file   : %10d\n', length(TIME));
fprintf('Number of layers          : %10d\n', NLAY);
fprintf('Number of Rows            : %10d\n', NROW);
fprintf('Number of columns         : %10d\n', NCOL);

if ~exist('userPeriods','var')
    userPeriods=periodsInFile;
else
    userPeriods=unique(userPeriods);
end

if any(userPeriods>max(periodsInFile))
    fprintf('You requested periods that are larger than maximum in file (%), last record will be used.\n',max(periodsInFile));
    userPeriods(userPeriods>max(periodsInFile))=max(periodsInFile);
    userPeriods=unique(userPeriods);
end

if userPeriods<1, C=[]; return; end

%% use user-specified transport steps or else all that are in the infile

if userPeriods==-1, userPeriods=max(periodsInFile); end

if  exist('userPeriods','var') && ~isempty(userPeriods)
    if any(userPeriods<min(periodsInFile)) || any(userPeriods>max(periodsInFile))
        error(['Requested periods beyond range(1..%d) in file %s!\n',...
               'Use userPeriods=-1 to select the last record\n',...
               'Use userPeriods< 1 to get the meta data of the file'],NTRP,fname);
    end
else
    userPeriods=unique(periodsInFile);
end

%% use user-specified layers or else all that are in the infile
if exist('userLayers'   ,'var') && ~isempty(userLayers)
    if ~isempty(userLayers(userLayers<1 | userLayers>NLAY))
        error('Requested layers are beyond range(1..%d) in file %s!',NLAY,fname);
    end
    userLayers=unique(userLayers);
else
    userLayers=unique(layersInFile);
end
LASTLAY=userLayers(end);

%% Selectors have values for all records of the input file.
%  Their contents is tranport-step index and the layer index.
%  The transport-step index is the number of the corresponding output record
%  The layer index that in the corresponding output 3D array
UTRPCOL=1; UPERCOL=2; ULAYCOL=3;

Select=zeros(nRec,ULAYCOL);
                             Select(:, UTRPCOL)   = trpStpsInFile;
for i=1:length(userPeriods), Select(periodsInFile == userPeriods(i) ,UPERCOL)=i; end % output rec index
for i=1:length(userLayers),  Select(layersInFile  == userLayers(i)  ,ULAYCOL)=i; end % ouput layer index

%% userLayers, userRows and userCols do not need an indicater array to pick out the
% correct data, the lists are used directly

if exist('userRows'   ,'var') && ~isempty(userRows)
    if ~isempty(userRows(userRows<1 | userRows>NROW))
        error('Requested userRows are beyond the range(1..%d) in file %s!',NROW,fname);
    end
else
    userRows=1:NROW;
end

if exist('userCols'   ,'var') && ~isempty(userCols)
    if ~isempty(userCols(userCols<1 | userCols>NCOL))
        error('Requested columns are beyond the range(1..%d) in file %s!',NCOL,fname);
    end
else
    userCols=1:NCOL;
end

%% What are the unique periods and tstep combinations? Because the number of
%  these equals the number of required output records

% Input records to deal with (unique transportstep-layer combinations)
IPTL=find(Select(:,UPERCOL)>0 & Select(:,ULAYCOL)>0 & Select(:,UTRPCOL)>0);

% Ouput records
nRecOut = length(find(Select(:,UPERCOL)>0 & Select(:,UTRPCOL)>0))/NLAY;

if rem(nRecOut,1)~=0,
    error('mflab:readMT3D:nRecOutNotAnInteger',...
        ['nRecOut (number of output rectorcs = %g, is not an integer !\n',...
        'Check if MT3D or SEAWAT finished normally.\n',...
        'Check your LST file, to see if something went wrong during the simulation.'],nRecOut);
end

B.values=NaN(length(userRows),length(userCols),length(userLayers)); % dimension of single outrec

C=repmat(B,[nRecOut,1]); % Size of output struct array is allocated here

NL=length(userLayers); iL=1; % counter and total number of user layers

%% Get the actual data values
iROut=1;
for i=1:length(IPTL)  % indices of required records in infile
    iRecIn=IPTL(i);
    fseek(fp,nbyt*(iRecIn-1),'bof'); % iROut starts at 0!

    iL   =Select(iRecIn,ULAYCOL);   % output layer  Nr

    [C(iROut).trpstp,...
        C(iROut).tstp,...
        C(iROut).period,...
        C(iROut).time,...
        C(iROut).label,...
        iLay,...
        C(iROut).rows,...
        C(iROut).cols,...
        values]=contrRec(fp,bytes);
    
    %fprintf('%4d %4d %4d %4d\n',iRecIn,iROut,iLay,iL)
    
    C(iROut).values(:,:,iL)=values(userRows,userCols);
    C(iROut).lays =userLayers;
    C(iROut).rows=userRows;
    C(iROut).cols=userCols;
    if verbose>0
        if iLay==LASTLAY
            fprintf('iRecIn=%4d trpStep=%4d, period=%3d, tstep=%3d, time=%12g, Layers=%3d --- iRecOut=%d\n',...
                iRecIn,...
                trpStpsInFile(iROut),...
                C(iROut).period,...
                C(iROut).tstp,...
                C(iROut).time,...
                length(C(iROut).lays),...
                iROut);
        end
    else
        if rem(i, 100)==0, fprintf('.'); end
        if rem(i,5000)==0, fprintf('iRecIn=%4d iRecOut=%4d\n',iRecIn,iROut); end
    end
    if iLay==LASTLAY, iROut=iROut+1; end
end
fprintf('%6d records read.\n.',iRecIn);
fprintf('%6d records in output struct.\n',iROut);
fprintf('\n');

fclose(fp);

end

function [ntrans,kstp,kper,time,label,ilay,nrow,ncol,values]=contrRec(fp,bytes)
% [ntrans,kstp,kper,time,label,ilay,nrow,ncol,data,readdataFlag]=contrRec(fid,bytes)
% --- reads a complete layer from a MODFLOW binary file
% --- if noread exists, only the info  record is read to safe time
% TO 070703

if nargin==1, % just get the byte offset in the Binary file
    fread(fp, 3,'int32');    % ntrans, kstp, kper
    fread(fp, 1,'float32');  % time
    
    % FORTRAN bytes at beginning and end of each record, are recognized by
    % assuming they are not in the file and checking wheather the first 4
    % bytes contain any non readable charachters if so, we assume that each
    % record is preceded and ended with 4 bytes
    % If compiled with a Fortran compiler that adds a different number of
    % bytes than 0 or 4 this has to be adapted
    % TO 100502
    
    p=ftell(fp);
    label =char(fread(fp, 4,'uchar')');
    if any(label<'0' & label~=' ') || any(label>'9' & label<'A') || any(label>'Z' & label<'a') || any(label>'z')
        bytes=length(label);
    else
        bytes=0;
    end
    fseek(fp,p,'bof');

    fread(fp,bytes,'int8');

    label =char(fread(fp,16,'uchar')');
    
    ncol  =fread(fp, 1,'int32');
    nrow  =fread(fp, 1,'int32');
    ilay  =fread(fp, 1,'int32');

    fread(fp,bytes,'int8');
    fread(fp,bytes,'int8');

    n=ftell(fp); fread(fp,1,'float'); floatlen=ftell(fp)-n;
    fseek(fp,floatlen*(ncol*nrow-1),0);
    
    fread(fp,bytes,'int8');

%    kstp=bytes;
    ntrans=bytes;
    return
end

% if nargin>1, then we know the number of bytes
    
fread(fp,bytes,'int8');

ntrans=fread(fp, 1,'int');
kstp  =fread(fp, 1,'int');
kper  =fread(fp, 1,'int');
time=fread(fp, 1,'float');

label =char(fread(fp,16,'char')');
label=label(label~=' ');

ncol  =fread(fp, 1,'int');
nrow  =fread(fp, 1,'int');
ilay  =fread(fp, 1,'int');

fread(fp,bytes,'int8');
fread(fp,bytes,'int8');

if nargout<9, % do not fetch the data
    n=ftell(fp); fread(fp,1,'float'); floatlen=ftell(fp)-n;
    fseek(fp,floatlen*(ncol*nrow-1),0);
else         % fetch the data
    values=permute(reshape(fread(fp,ncol*nrow,'float'),[ncol,nrow]),[2,1]);
end

fread(fp,bytes,'int8');

end
