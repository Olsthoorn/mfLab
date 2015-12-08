function HOut=readDat(fName,userPeriods,userTSteps,userLays,userRows,userCols)
%READDAT reads binary MODFLOW heads or drawdown output into a struct, H
%
% Example:
%   H=readDat([+]fName [,userPeriods [,tsteps [,userLays [,userRows [,userCols]]]]]])
%
%   H contains both the meta data and the actual values.
%   Use the optional arguments to specify your selection or omit the options
%
%   to get all (by just specifying the filename).
%   Prepend filename with a "+" to get verbose screen output.
%   e.g.: H=readat(['+' fName '.HDS']);
%
%   To get only the meta data of all headers (the meta data) and printout
%   the overal meta data of the file do this:
%   H=readDat(fName,-1)
%
%   To get all data from file do this
%   H=readDat(fName);
%
%   Use the optional arguments to make subselections from this file
%   which will reduce the required workspace memory and time to read the file
%   Optional arguments may be omittted as shown or used empty '' or [].
%   So, to fetch all stress periods, but only time setps 1,3 and 6,all layers,
%   only row 3 and all columns, do this:
%   H=readDat(fName,'',[1,3,6],'',3)
%   (you don't need the last argument if left open.

% To fetch stress periods 1,3 and 5..8, all time steps, layer 3,
% rows 3..10 and columns 12..21, do this:%
%   H=readDat(fName,[1 3 5:8],'',5,3:10,12:21}
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
%    .time         same, for convenience
%    .NCOL         number of original columns in file fName
%    .NROW         number of original rows    in file fName
%    .NLAY         number of original layers  in file fName
%    .Cols         list of fetched columns
%    .Rows         list of fetched rows
%    .Lays         list of fetched layers
%
% However, if just the metadata is requrested using -1 as period option,
% the struct H contains headers of all records in the file. Its fields are
%    .tstp       time step
%    .period         stress period
%    .pertim      time elapsed in current time stress period
%    .totim       time elapsed since start of simulation
%    .time        same as totim
%    .label       label of this record
%    .iLay        layer number of this record
%    .NROW        number of rows
%    .NCOL        number of columns
%
%  SEE ALSO: readBud, readMT3D
%
%   TO 090104 091214
%   TO 100818 logic redone to clean up necessary to fix number of ouput recs
%   TO 130405 added time and changed tstep to tstp for conpatibility with readBud and readMT3D
%
% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

% mask value, so set inactive cells to NaN including HDry cells
maskVal = 1e9;

if     fName(1)=='+' || fName(1)=='1', fName=fName(2:end); verbose=1;
elseif fName(1)=='-' || fName(1)=='0', fName=fName(2:end); verbose=0;
else   verbose=0;
end

doAllPeriods = ~exist('userPeriods','var')  || isempty(userPeriods);
doAllTstps   = ~exist('userTSteps'  ,'var') || isempty(userTSteps);
doAllLayers  = ~exist('userLays' ,'var') || isempty(userLays);
doAll        = doAllPeriods && doAllTstps && doAllLayers;

fp=fopen(fName); if fp<1, error('READDAT: Can''t find or open file %s!',fName); end

fprintf('Reading MODFLOW binary output file <<%s>> verbose= %d\n',fName,verbose);

fseek(fp,0,1); n=ftell(fp); fseek(fp,0,-1);
if n==0, error('mflab:readDat:filehaszeorlength',...
        ['%s: file %s has zero length, this normally means that the model has not\n',...
         'finished normally.\n',...
         'REMEDY: Check the LST file to see what went wrong.'],mfilename,fName);
end

%% Advance a single cell-by-cell record & compute the number or records in file

bytes=contrRec(fp);

nbyt=ftell(fp); fseek(fp,0,'eof'); fLen=ftell(fp); nRec=fLen/nbyt; % number of records in file

if rem(fLen,nbyt)~=0,
    error('BINARY file %s not of standard type. Sorry, can''t read it, use a better compiler!',fName);
end

%% Get the information about the contents of the file from record headings

fprintf('Scanning headers\n');

H=repmat(...
    struct(...    
        'tstp',NaN,...
        'period',NaN,...
        'pertim',NaN,...
        'totim',NaN,...
        'time',NaN,...
        'label',{[]},...
        'iLay',NaN,...
        'NROW',NaN,...
        'NCOL',NaN,...
        'iROut',NaN,...
        'mustSelect',0),...
    nRec,1);

%% get meta data for each record in file

nRecOut=0; oldPeriod=-1; oldTStp=-1;

for i=1:nRec
    fseek(fp,(i-1)*nbyt,'bof');  % move to start of last layer in file
    
    [H(i).tstp,...
        H(i).period,...
        H(i).pertim,...
        H(i).totim,...
        H(i).label,...
        H(i).iLay,...
        H(i).NROW,...
        H(i).NCOL...
        ]=contrRec(fp,bytes,fName);
    
        %% Do we need this record, if so, set output record number
    % We have foreach(SP) --> timesteps --> foreach(timestep) --> labels
    
    newPeriod= H(i).period  ~=oldPeriod;
    newTStp  = H(i).tstp    ~=oldTStp;
    
    H(i).mustSelect= doAll;
    if ~doAll
        H(i).mustSelect    = doAllPeriods || any(H(i).period==userPeriods);
        if H(i).mustSelect
           H(i).mustSelect = doAllTstps   || any(H(i).tstp  ==userTSteps);
           if H(i).mustSelect
               H(i).mustSelect = doAllLayers || any(H(i).iLay==userLays);
           end
        end
    end
      
    if H(i).mustSelect,
        if newPeriod || newTStp,
            nRecOut=nRecOut+1;
        end
        H(i).iROut=nRecOut;
        oldTStp   = H(i).tstp;
        oldPeriod = H(i).period;
    end
    
    if rem(i, 100)==0, fprintf('.'); end
    if rem(i,5000)==0, fprintf('%d records read\n',i); end
end

NLAY=max([H.iLay]);

if ~exist('userLays','var') || isempty(userLays), userLays=1:NLAY;        end
if ~exist('userRows','var') || isempty(userRows), userRows=1:H(end).NROW; end
if ~exist('userCols','var') || isempty(userCols), userCols=1:H(end).NCOL; end

fprintf('finished, %d records scanned\n',i);

if exist('userPeriods','var') && ~isempty(userPeriods) && userPeriods(1)==0
    return;  % just return the header data
end

%% Evaluate what we've got in terms of the userPeriods, time steps unique lables etcetera before we continue to select the actual data

NPER=max([H.period]);
NSTP=max([H.tstp]);
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


%% Data struct

ISelect= find([H.mustSelect]);

LASTLAY=userLays(end);

HOut=repmat(...
    struct(...
       'tstp',NaN,...
        'period',NaN,...
        'pertim',NaN,...
        'totim', NaN,...
        'time',NaN,...
        'label',{[]},...
        'NROW',NaN,...
        'NCOL',NaN,...
        'values',NaN(length(userRows),length(userCols),length(userLays)),...
        'lays',userLays,...
        'rows',userRows,...
        'cols',userCols),...
     nRecOut,1);

%% Get the actual data values

fprintf('Reading requested data ...\n');

iL=0; oldPeriod=-1; oldTStp=-1;
for i=1:length(ISelect)  % indices of required records in infile
    iRecIn=ISelect(i);
    
    fseek(fp,nbyt*(iRecIn-1),'bof');
    
    iROut = H(iRecIn).iROut;
    
   [HOut(iROut).tstp,...
    HOut(iROut).period,...
    HOut(iROut).pertim,...
    HOut(iROut).totim,...
    HOut(iROut).label,...
    iLay,...         % the layer number of this record
    HOut(iROut).NROW,...
    HOut(iROut).NCOL,...
    values...        % entire layer array
    ]=contrRec(fp,bytes,fName);

    %fprintf('%4d %4d %4d %4d\n',iRecIn,iROut,iLay,iL)
        
    newPeriod = H(iRecIn).period~= oldPeriod;
    newTStp   = H(iRecIn).tstp ~= oldTStp;
    
    if newPeriod|| newTStp,
        iL=1;
    else
        iL=iL+1;
    end

    oldPeriod = H(iRecIn).period;
    oldTStp   = H(iRecIn).tstp;
    
    HOut(iROut).values(:,:,iL)=values(userRows,userCols);
    HOut(iROut).lays=userLays;
    HOut(iROut).rows=userRows;
    HOut(iROut).cols=userCols;
    if verbose>0
        if iLay==LASTLAY,  % every time iL corresponds to the last user layer
            fprintf('iRecIn=%4d period=%3d, tstp=%3d, pertim=%12g, totim=%12g Layers=%3d iRecOut=%4d\n',...
                iRecIn,...
                HOut(iROut).period,...
                HOut(iROut).tstp,...
                HOut(iROut).pertim,...
                HOut(iROut).totim,...
                length(HOut(iROut).lays),...
                iROut);
        end
    else
        if rem(i, 100)==0, fprintf('.'); end
        if rem(i,5000)==0, fprintf('iRecIn=%d iRecOut=%d\n',iRecIn,iROut); end
    end
    
    %% Mask HOut.values
    % TO 130425
    
    HOut(iROut).values(HOut(iROut).values>maskVal)=NaN;
    
%     if iL==LASTLAY
%         iROut=iROut+1;
%     end
end
fprintf('%6d records read.\n.',iRecIn);
fprintf('%6d records in output struct.\n',iROut);
fprintf('\n');

fclose(fp);

%% add time for convenience, so that same as in Conc and Bud  TO 130405
for it = 1:numel(HOut)
    HOut(it).time = HOut(it).totim;
end



end

function [kstp,kper,pertim,totim,label,nlay,nrow,ncol,values]=contrRec(fp,bytes,fName)
% [kstp,kper,text,nlay,nrow,ncol,data]=contrRec(fp,bytes[,fName])
% --- reads a complete layer from a MODFLOW binary file
% --- if noread exists, only the info  record is read to safe time
% TO 070703 091214

n = ftell(fp); fseek(fp,0,1); fLen = ftell(fp); fseek(fp,n,-1);

if nargin==1, % just get the byte offset in the Binary file
    
    kstp  =fread(fp, 1,'int32'); %#ok
    kper  =fread(fp, 1,'int32');

    pertim=fread(fp, 1,'float32');
    totim =fread(fp, 1,'float32');
    
    % FORTRAN bytes at beginning and end of each record, are recognized by
    % assuming they are not in the file and checking wheather the first 4
    % bytes contain any non readable charachters if so, we assume that each
    % record is preceded and ended with 4 bytes
    % If compiled with a Fortran compiler that adds a different number of
    % bytes than 0 or 4 this has to be adapted
    % TO 100502
    
    p=ftell(fp);
%    label =char(fread(fp, 4,'uchar')');
    label =char(fread(fp,20,'uchar')');
%    bytes = regexp(label,' +[A-Z0-9]+')-1;
    bytes = regexp(label,'    +|[A-Z0-9]{4,}')-1;  % TO 120925
    bytes = bytes(1);

    fseek(fp,p,'bof');

    try
        fread(fp,bytes,'int8');
    catch ME
        fprintf(2,['Can''t read to required end of HDS or DDN file\n',...
            'Check your LST file to see that your MODFLOW has terminated normally.\n',...
            'See that it did not fail to converge, so that it did couldn''t complete the file !\n']);
        throw(ME);
    end
    
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
    if ftell(fp)>=fLen
        error('trying to read beyond end of file %s, fLen = %d',fName,fLen);
    end
    values=permute(reshape(fread(fp,ncol*nrow,'float'),[ncol,nrow]),[2,1]);    
end

fread(fp,bytes,'int8');

end

