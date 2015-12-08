function COut=readMT3D(varargin)
%READMT3D reads unformatted MT3DMS concentration file into struct 
%
% Example:
%   C=readMT3D([-]fname [,userPeriods [,userLays [,userRows [,userCols]]]]]][,options])
%   C=readMT3D(fname,-1); % to get info on the contents of the file.
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
% options: this is recognized
%   C=readMT3D(... 'round',digits,'min',minVal,'max',maxVal ...);
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
% TO 120205 logic redone to speed up selection from large data files
%
% Copyright 2009-2012 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later
[digits,varargin] = getProp(varargin,'round',[]);
[maxVal,varargin] = getProp(varargin,'max'  ,[]);
[minVal,varargin] = getProp(varargin,'min'  ,[]);

[fname, varargin] = getNext(varargin,'char','MT3D000.UCN');

if     fname(1)=='+' || fname(1)=='1', fname=fname(2:end); verbose=1;
elseif fname(1)=='-' || fname(1)=='0', fname=fname(2:end); verbose=0;
else   verbose=0;
end

[userPeriods, varargin] = getNext(varargin,'double',[]);
[userTSteps,  varargin] = getNext(varargin,'double',[]);
[userLays,    varargin] = getNext(varargin,'double',[]);
[userRows,    varargin] = getNext(varargin,'double',[]);
[userCols,    varargin] = getNext(varargin,'double',[]);

if ~isempty(varargin)
    msgId = 'mfLab:readMT3D:unusedArguments';
    warning('on',msgId)
    warning('Some arguments have not been used in the call to %s',mfilename);
    warning('off',msgId);
end

doAllPeriods = isempty(userPeriods);
doAllTSteps  = isempty(userTSteps);
doAllLayers  = isempty(userLays);
doAll        = isempty(doAllPeriods) && isempty(doAllLayers);

fp=fopen(fname); if fp<1, error('%s: Can''t find or open file %s!',mfilename,fname); end

fprintf('Reading MT3DMS binary output file <<%s>> verbose= %d\n',fname,verbose);

%% Advance a single cell-by-cell record & compute the number or records in file

bytes=contrRec(fp);  % bytes is compiler specific UNFORMATTED file record heading

nbyt=ftell(fp); fseek(fp,0,'eof'); Nbyt=ftell(fp); nRec=Nbyt/nbyt; % number of records in file

if rem(Nbyt,nbyt),
    warning('mfLab:readMT3D:fileLengthDoesNotMatchREcordLength',...
        ['The UNFORMATTED MT3D output file %s seems to be of non standard type.\n'...
        'or incomplete. Its total length does not match an integer number of layer records.\n',...
        'The file may be corrupt or incomplete, maybe because MT3DMS or SEAWAT stopped prematurely.\n',...
        'Less probably it might also be due to the fortran code not producing binary output files\n',...
        'of standard. In that case you may need to recompile the source code with another FORTRAN compiler.\n',...
        'I continue as far a possible. If the majority of the output is correct (see last records),\n',...
        ' then the reason will be the first, premature stop of the groundwater code.]'],fname,nRec);
end

nRec = floor(nRec);

%% Get the information about the contents of the file from record headings

fprintf('Scanning headers...\n');

C=repmat(...
    struct(...    
        'trpstp',NaN,...
        'period',NaN,...
        'tstp',NaN,...
        'time',NaN,...
        'label',{[]},...
        'NLAY',NaN,...
        'NROW',NaN,...
        'NCOL',NaN,...
        'iROut',NaN,...
        'mustSelect',0),...
    nRec,1);

%% Get headers

nRecOut=0; oldPeriod=-1; oldTStp=-1; oldTrpStp=-1;

for i=1:nRec
    fseek(fp,(i-1)*nbyt,'bof');  % move to start of last layer in file
    
    [C(i).trpstp,...
        C(i).tstp,...
        C(i).period,...
        C(i).time,...
        C(i).label,...
        C(i).iLay,...
        C(i).NROW,...
        C(i).NCOL...
        ]=contrRec(fp,bytes);

    newPeriod= C(i).period  ~=oldPeriod;
    newTStp  = C(i).tstp    ~=oldTStp;
    newTrpStp= C(i).trpstp  ~=oldTrpStp;
    
    C(i).mustSelect= doAll;
    if ~doAll
        C(i).mustSelect    = doAllPeriods || any(C(i).period==userPeriods);
        if C(i).mustSelect
           C(i).mustSelect = doAllTSteps   || any(C(i).tstep  ==userTSteps);
           if C(i).mustSelect
               C(i).mustSelect = doAllLayers || any(C(i).iLay==userLays);
           end
        end
    end
      
    if C(i).mustSelect,
        if newPeriod || newTStp || newTrpStp,
            nRecOut=nRecOut+1;
        end
        C(i).iROut=nRecOut;
        oldTStp   = C(i).tstp;
        oldPeriod = C(i).period;
        oldTrpStp = C(i).trpstp;
    end

    if rem(i, 100)==0, fprintf('.'); end
    if rem(i,5000)==0, fprintf('%d records read\n',i); end

end
fprintf('...finished, %d records read\n',i);

%% Evaluate what we've got in terms of the periods, timeInFile steps unique lables etcetera before we continue to select the actual data

NPER = max([C.period]);
NSTP = max([C.tstp]);
NLAY = max([C.iLay]);
NROW = C(end).NROW;
NCOL = C(end).NCOL;
TIME = unique([C.time]);

if isempty(userLays),    userLays=1:NLAY; end
if isempty(userRows),    userRows=1:NROW; end
if isempty(userCols),    userCols=1:NCOL; end
if isempty(userPeriods), userPeriods=1:NPER; end
if isempty(userTSteps),  userTSteps=1:NSTP; end

fprintf('File contains the following:\n');
fprintf('Number of records in file : %10d\n', nRec);
fprintf('Number of stress periods  : %10d\n', NPER);
fprintf('Number of time steps      : %10d\n', NSTP);
fprintf('Highest time in file      : %10g\n', max(TIME));
fprintf('Number of times in file   : %10d\n', length(TIME));
fprintf('Number of layers          : %10d\n', NLAY);
fprintf('Number of Rows            : %10d\n', NROW);
fprintf('Number of columns         : %10d\n', NCOL);

if any(userPeriods>NPER)
    error('userPeriods must be <= NPER=%d',NPER);
end
if any(userTSteps>NSTP)
    error('userTSteps must be <= NSTP=%d\n',NSTP);
end

if userPeriods<1; return; end

if rem(nRecOut,1)~=0,
    error('mflab:readMT3D:nRecOutNotAnInteger',...
        ['nRecOut (number of output rectorcs = %g, is not an integer !\n',...
        'Check if MT3D or SEAWAT finished normally.\n',...
        'Check your LST file, to see if something went wrong during the simulation.'],nRecOut);
end


%% Struct to store actual data

COut=repmat(...
    struct(...
        'trpstp',NaN,...
        'tstp',  NaN,...
        'period',NaN,...
        'time',  NaN,...
        'label', {[]},...
        'lays',  NaN,...
        'rows',  NaN,...
        'cols',  NaN,...
        'values',NaN(length(userRows),length(userCols),length(userLays))...
    ), nRecOut,1);

%% Get the actual data values
ISelect = find([C.mustSelect]);

LASTLAY=userLays(end);

iL=0; oldPeriod=-1; oldTStp=-1;

for i=1:length(ISelect)  % indices of required records in infile
    iRecIn=ISelect(i);
    
    fseek(fp,nbyt*(iRecIn-1),'bof'); % iRecIn starts at 0!

    iROut=C(iRecIn).iROut;

    [COut(iROut).trpstp,...
        COut(iROut).tstp,...
        COut(iROut).period,...
        COut(iROut).time,...
        COut(iROut).label,...
        iLay,...
        COut(iROut).rows,...
        COut(iROut).cols,...
        values]=contrRec(fp,bytes);
    
    %fprintf('%4d %4d %4d %4d\n',iRecIn,iROut,iLay,iL)
    
    newPeriod = C(iRecIn).period ~= oldPeriod;
    newTStp   = C(iRecIn).tstp   ~= oldTStp;
    newTrpStp = C(iRecIn).trpstp ~= oldTrpStp;

    if newPeriod|| newTStp || newTrpStp,
        iL=1;
    else
        iL=iL+1;
    end

    oldPeriod = C(iRecIn).period;
    oldTStp   = C(iRecIn).tstp;
    oldTrpStp = C(iRecIn).trpstp;
    
    COut(iROut).values(:,:,iL)=values(userRows,userCols);
    COut(iROut).lays =userLays;
    COut(iROut).rows=userRows;
    COut(iROut).cols=userCols;
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
%     if iLay==LASTLAY, iROut=iROut+1; end
end

%% Masking options
if ~isempty(digits)
    for it=1:numel(COut)
        COut(it).values = round(COut(it).values,digits);
    end
end

if ~isempty(minVal)
    for it=1:numel(COut)
        COut(it).values = max(COut(it).values,minVal);
    end
end

if ~isempty(maxVal)
    for it=1:numel(COut)
        COut(it).values = min(COut(it).values,maxVal);
    end
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
%    label =char(fread(fp, 4,'uchar')');
    label =char(fread(fp,20,'uchar')');
    bytes = regexp(label,'CONCENTRATION')-1;
%     if any(label<'0' & label~=' ') || any(label>'9' & label<'A') || any(label>'Z' & label<'a') || any(label>'z')
%         bytes=length(label);
%     else
%         bytes=0;
%     end
    fseek(fp,p,'bof');

    try
        fread(fp,bytes,'int8');
    catch ME
        fprintf(2,['Can''t read to required end of MT3D concentration file\n',...
            'Check your LST file to see that MT3D/SEAWAT has terminated normally.\n',...
            'See that it did not fail to converge, so that it did couldn''t complete the file !\n']);
        throw(ME);
    end
    
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
