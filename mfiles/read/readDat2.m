function H=readDat2(fname,userPeriods,userTstp,userLayers,userRows,userCols)
%READDAT2 reads binary MODFLOW heads or drawdown output into a struct, H
%
% USAGE:
%   H=readDat([+]fname [,userPeriods [,tsteps [,userLayers [,userRows [,userCols]]]]]])
%
%   H contains both the meta data and the actual values.
%   Use the optional arguments to specify your selection or omit the options
%
% EXAMPLE:
%   to get all (by just specifying the filename).
%   Prepend filename with a "+" to get verbose screen output.
%   e.g.: H=readat(['+' fname '.HDS']);
%
%   To get only the meta data of all headers (the meta data) and printout
%   the overal meta data of the file do this:
%   H=readDat(fname,-1)
%
%   To get all data from file do this
%   H=readDat(fname);
%
%   Use the optional arguments to make subselections from this file
%   which will reduce the required workspace memory and time to read the file
%   Optional arguments may be omittted as shown or used empty '' or [].
%   So, to fetch all stress periods, but only time setps 1,3 and 6,all layers,
%   only row 3 and all columns, do this:
%   H=readDat(fname,'',[1,3,6],'',3)
%   (you don't need the last argument if left open.

% To fetch stress periods 1,3 and 5..8, all time steps, layer 3,
% rows 3..10 and columns 12..21, do this:%
%   H=readDat(fname,[1 3 5:8],'',5,3:10,12:21}
%
% OUTPUT:
%   H is a cell array whose length is the number of records in the file
%   it has the following fields:
%    .values is a 3D array containing the requested concentration values
%    .label        label telling the contents
%    .period       stress period number
%    .tstp         time step number within stress period
%    .pertim       time since begin of stress period
%    .totim        time since begin of simulation
%    .NCOL         number of original columns in file fname
%    .NROW         number of original rows    in file fname
%    .NLAY         number of original layers  in file fname
%    .Cols         list of fetched columns
%    .Rows         list of fetched rows
%    .Lays         list of fetched layers
%
% However, if just the metadata is requrested using -1 as period option,
% the struct H contains headers of all records in the file. Its fields are
%    .tstep       time step
%    .per         stress period
%    .pertim      time elapsed in current time stress period
%    .totim       time elapsed since start of simulation
%    .label       label of this record
%    .iLay        layer number of this record
%    .NROW        number of rows
%    .NCOL        number of columns
%
%  SEE ALSO: readBud, readMT3D
%
%   TO 090104 091214
%   TO 100818 logic redone to clean up necessary to fix number of ouput recs
%   TO 120205 logic completely redone to speed up selection from large files
%
% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if     fname(1)=='+' || fname(1)=='1', fname=fname(2:end); verbose=1;
elseif fname(1)=='-' || fname(1)=='0', fname=fname(2:end); verbose=0;
else   verbose=0;
end

fp=fopen(fname); if fp<1, error('READDAT: Can''t find or open file %s!',fname); end

fprintf('Reading MODFLOW binary output file <<%s>> verbose= %d\n',fname,verbose);

fseek(fp,0,1); n=ftell(fp); fseek(fp,0,-1);
if n==0, error('mflab:readDat:filehaszeorlength','file %s has zero length',fname); end

%% Advance a single cell-by-cell record & compute the number or records in file

bytes=contrRec(fp);

nbyt=ftell(fp); fseek(fp,0,'eof'); Nbyt=ftell(fp); nRec=Nbyt/nbyt; % number of records in file

if rem(Nbyt,nbyt)~=0,
    error('BINARY file %s not of standard type. Sorry, can''t read it, use a better compiler!',fname);
end

%% Get the information about the contents of the file from record headings

fprintf('Scanning headers\n');
H=[];
H(nRec).tstep=NaN;
H(nRec).per=NaN;
H(nRec).pertim=NaN;
H(nRec).totim=NaN;
H(nRec).label='NUL';
H(nRec).iLay=NaN;
H(nRec).NROW=NaN;
H(nRec).NCOL=NaN;
for i=1:nRec  % get meta data for each record in file
    fseek(fp,(i-1)*nbyt,'bof');  % move to start of last layer in file
    
    [H(i).tstep,...
        H(i).per,...
        H(i).pertim,...
        H(i).totim,...
        H(i).label,...
        H(i).iLay,...
        H(i).NROW,...
        H(i).NCOL...
        ]=contrRec(fp,bytes);
    
    if rem(i, 100)==0, fprintf('.'); end
    if rem(i,5000)==0, fprintf('%d records read\n',i); end
end

fprintf('finished, %d records read\n',i);
if exist('userPeriods','var') && ~isempty(userPeriods) && userPeriods(1)==0
    return;  % just return the header data
end

%% Evaluate what we've got in terms of the userPeriods, time steps unique lables etcetera before we continue to select the actual data

NPER=max([H.per]);
NSTP=max([H.tstep]);
NLAY=max([H.iLay]);
NROW=max([H.NROW]);
NCOL=max([H.NCOL]);
totim=max([H.totim]);

fprintf('File contains the following:\n');
fprintf('Number of records in file  : %10d\n',nRec);
fprintf('Number of stress userPeriods: %10d\n',NPER);
fprintf('Number of time steps       : %10d\n',NSTP);
fprintf('Number of layers           : %10d\n',NLAY);
fprintf('Number of rows             : %10d\n',NROW);
fprintf('Number of columns          : %10d\n',NCOL);
fprintf('Maximum time in file       " %10g\n',totim);

if exist('userPeriods','var') && ~isempty(userPeriods) && userPeriods(1)==0,
    return;  % just return the meta data
end

%% select the speciied zones or all if nothing specified
periodsInFile=[H.per]';   %periodsInFile  the period numbers in the InFile
tstepsInFile =[H.tstep]'; %tstepsInFile   the time step numbers in the InFile
layersInFile =[H.iLay]';  %layersInFile   the layer numbefs in the InFile

%% the period numbers requested by the user or all of them
if ~exist('userPeriods','var'),
    userPeriods=periodsInFile;
else
    userPeriods=unique(userPeriods);  % this is set of user requrested stress userPeriods
end

if any(userPeriods>max(periodsInFile))
    fprintf('You requrested periods >max(periodsInFile) (=%d), last record will be used\n',max(periodsInFile));
    userPeriods(userPeriods>max(periodsInFile))=max(periodsInFile);
end

if  any(userPeriods<min(periodsInFile)) || any(userPeriods>max(periodsInFile))
    error(['Requested userPeriods beyond range(1..%d) in dat file %s!\n'...
           'Use userPeriods=-1 to only the last record\n',...
           'Use userPeriods= 0 to get the metaData of the file.'],NPER,fname);
end

%% the time step numbers requested by the user or all of them
if exist('userTstp'  ,'var') && ~isempty(userTstp)  % if the user specified which time steps he wants
    if ~isempty(userTstp(userTstp<1 | userTstp>NSTP))
        error('Requested userPeriods beyond range(1..%d) in buget file %s!',NSTP,fname);
    end
    userTstep=unique(userTstp);
else
    userTstep=unique(tstepsInFile);
end

%% the layer numbers requested by the user of all of them
if exist('userLayers'   ,'var') && ~isempty(userLayers)  % if the user specifies which layers he/she wants
    if ~isempty(userLayers(userLayers<1 | userLayers>NLAY))
        error('Requested layers are beyond range(1..%d) in dat file %s!',NLAY,fname);
    end
    userLayers=unique(userLayers);
else
    userLayers=unique(layersInFile);
end
LASTLAY=userLayers(end);

%% Selectors have values for all records of the inptut file. There
%  contents are userPeriods index, the userTstep index, the userLayers index
%  the combined userPeriods and userTstep index plus
%  the output record index and the output layer index
UPCOL=1; UTCOL=2; ULCOL=3; UPTCOL=4; OLCOL=6;

fprintf('...please wait, I''m busy setting up selectors ...\n');

Select=zeros(nRec,OLCOL);  
for i=1:length(userPeriods), Select([H.per]  ==userPeriods(i),UPCOL)=i; end
for i=1:length(userTstep),  Select([H.tstep]==userTstep(i) ,UTCOL)=i; end
for i=1:length(userLayers), Select([H.iLay] ==userLayers(i),ULCOL)=i; end

% input records to deal with
IPT=find(Select(:,UPCOL)>0 & Select(:,UTCOL)>0 & Select(:,ULCOL));

% output records
userPT = unique([periodsInFile(IPT) tstepsInFile(IPT)],'rows');

for i=1:size(userPT,1), Select(periodsInFile ==userPT(i,1) & tstepsInFile==userPT(i,2),UPTCOL)=i; end

%% rows and cols do not need an indicater array to pick out the correct
%  data, the lists are used directly

if exist('userRows'   ,'var') && ~isempty(userRows)
    if ~isempty(userRows(userRows<1 | userRows>NROW))
        error('Requested rows are beyond the range(1..%d) in dat file %s!',NROW,fname);
    end
else
    userRows=1:NROW;
end

if exist('userCols'   ,'var') && ~isempty(userCols)
    if ~isempty(userCols(userCols>NCOL | userCols<1))
        error('Requested columns are beyond the range(1..%d) in dat file %s!',NCOL,fname);
    end
else
    userCols=1:NCOL;
end

%% What are the unique periods and tstep combinations? Because the number of

nRecOut = size(userPT,1);

if rem(nRecOut,1)~=0,
    error('mflab:readDat:nRecOutNotAnInteger',...
        ['nRecOut (number of output rectorcs = %g, is not an integer !\n',...
        'Check if MT3D or SEAWAT finished normally.\n',...
        'Check your LST file, to see if something went wrong during the simulation.'],nRecOut);
end


B.values=NaN(length(userRows),length(userCols),length(userLayers)); % dimension of single outrec

H=repmat(B,[nRecOut,1]); % size of output sttruct array is allocated here

%% Get the actual data values

iROut=1;
fprintf('Reading requested data ...\n');
for i=1:length(IPT)  % indices of required records in infile
    iRecIn=IPT(i);
    fseek(fp,nbyt*(iRecIn-1),'bof');
    %iROut=Select(iRecIn,UPTCOL);   % output record Nr
    iL=Select(iRecIn,ULCOL);    % output layer  Nr
    
   [H(iROut).tstp,...
    H(iROut).period,...
    H(iROut).pertim,...
    H(iROut).totim,...
    H(iROut).label,...
    iLay,...         % the layer number of this record
    H(iROut).NROW,...
    H(iROut).NCOL,...
    values...        % entire layer array
    ]=contrRec(fp,bytes);

    %fprintf('%4d %4d %4d %4d\n',iRecIn,iROut,iLay,iL)
        
    H(iROut).values(:,:,iL)=values(userRows,userCols);
    H(iROut).lays=userLayers;
    H(iROut).rows=userRows;
    H(iROut).cols=userCols;
    if verbose>0
        if iLay==LASTLAY,  % everty time iL corresponds to the last user layer
            fprintf('iRecIn=%4d per=%3d, tstep=%3d, pertim=%12g, totim=%12g Layers=%3d iRecOut=%4d\n',...
                iRecIn,...
                H(iROut).period,...
                H(iROut).tstp,...
                H(iROut).pertim,...
                H(iROut).totim,...
                length(H(iROut).lays),...
                iROut);
        end
    else
        if rem(i, 100)==0, fprintf('.'); end
        if rem(i,5000)==0, fprintf('iRecIn=%d iRecOut=%d\n',iRecIn,iROut); end
    end
    if iL==LASTLAY
        iROut=iROut+1;
    end
end
fprintf('%6d records read.\n.',iRecIn);
fprintf('%6d records in output struct.\n',iROut);
fprintf('\n');

fclose(fp);

end

function [kstp,kper,pertim,totim,label,nlay,nrow,ncol,values]=contrRec(fp,bytes)
% [kstp,kper,text,nlay,nrow,ncol,data]=contrRec(fp,bytes)
% --- reads a complete layer from a MODFLOW binary file
% --- if noread exists, only the info  record is read to safe time
% TO 070703 091214

if nargin==1, % just get the byte offset in the Binary file
    fread(fp, 2,'int32');     % kstp, kper
    fread(fp, 2,'float32');   % pertim, totim
    
    
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
    nlay  =fread(fp, 1,'int32');

    fread(fp,bytes,'int8');
    fread(fp,bytes,'int8');

    n=ftell(fp); fread(fp,1,'float'); floatlen=ftell(fp)-n;
    fseek(fp,floatlen*(ncol*nrow-1),0);
    
    fread(fp,bytes,'int8');

    kstp=bytes;
    return
end


% if bytes are known, nargin>1

fread(fp,bytes,'int8');

kstp  =fread(fp, 1,'int32');
kper  =fread(fp, 1,'int32');

pertim=fread(fp, 1,'float32');
totim =fread(fp, 1,'float32');

label =char(fread(fp,16,'uchar')'); label=label(label~=' ');

ncol  =fread(fp, 1,'int32');
nrow  =fread(fp, 1,'int32');
nlay  =fread(fp, 1,'int32');

fread(fp,bytes,'int8');
fread(fp,bytes,'int8');

if nargout<9 % don't fetch
    n=ftell(fp); fread(fp,1,'float32'); floatlen=ftell(fp)-n;
    fseek(fp,floatlen*(ncol*nrow-1),0);
else
    values=permute(reshape(fread(fp,ncol*nrow,'float'),[ncol,nrow]),[2,1]);    
end

fread(fp,bytes,'int8');

end

