function [B,budlabels]=readBud2(fname,userLabels,userPeriods,userTstps,userLayers,userRows,userCols)
%READBUD2 reads MODFLOW budget file into a Matlab structure array
%
% superseded by readBud
%
% USAGE:
%   [B,budlabels]=readBud([+]fname [,userLabels [,userPeriods [,tsteps [,userLayers [,userRows [,userCols]]]]]])
%
%   Fast select any portion of the budget file, no matter how big the file is.
%
%   budlabels are the list of unique labels in the entire budget file
%
% Example1:
%    To get only the headers, i.e. meta data from the file, do this:
%    B=readBud(fname,-1);
%
% Examle2:
%    to fetch all data from file do this
%    B=readBud(fname);
%
% Example3:
%    to get more verbose screen output prepend fname with a "+" 
%    B=readBud(['+' fname]);
%
% Example3:
%    Use the optional arguments to subselect data from file. For instance:
%    extracts all userLabels, but only stress userPeriods 1,3 and 6 and only
%    layer 3:
%    B=readbud(fname,'',[1,3,6],'',3)
%
% Example5:
%    To extract only the 'STORAGE' and 'FLOWRIGHTFACE' flow terms and these only
%    for stress periods 1,3 and 6, all time steps and layer 5 and rows 3...10
%    and columns 12...31, do this:
%    readbud(fname,{'STORAGE','FLOWRIGHTFACE'},[1 3 5],'',5,3:10,12:21),
%
% If only 1 flow term is desired you may omit the braces { }.
% The layers, rows and columns can be extracted in any order by adequate
% standard Matlab indexing.
%
% The resulting variable, here named B, is an struct array. Its length
% equals the number of records in the file (or leff if you limited your
% selection by qualifying userPeriods and or timesteps).
% This struct has the following fields
%   .label        cell array of the label names in the file, which refer to
%                 the cell by cell flow terms it contains
%   .term(NLabel) cell array with same length and order as label in which
%                 each cell has a full 3D array cell by cell flow terms
%                 corresponding to the label in the same position
%   .pertim       time within this period
%   .totim        total time from start of simulation
%   .period       current stress period number
%   .tstp         current time step number in this stress period
%   .NROW         total number of userRows    in budget file
%   .NCOL         total number of columns in budget file
%   .NLAY         total number of layers  in budget file
%   .cols         list of the extracted columns
%   .rows         list of the extracted userRows
%   .lays         list of the extracted layers
%   
% To get to a row and column in a layer of a given flow term
% that corresponds to the second label, do this:
% B(iRecord).term{iLabel}(yourrows,yourcols,yourlays)
% which is the standard Matlab way of indexing.
%
% SEE ALSO: readBud readDat readMT3D
%
% TO 090104 091214
% TO 100602 make budlabels unique because labels may vary between stress userPeriods
% TO 100818 logic redone to clean up in accordance with readDat and readMT3D
% TO 100917 added bulabels to output
%

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later


if     fname(1)=='+' || fname(1)=='1', fname=fname(2:end); verbose=1;
elseif fname(1)=='-' || fname(1)=='0', fname=fname(2:end); verbose=0;
else   verbose=0;
end

fp=fopen(fname); if fp<1, error('READBUD: Can''t find or open file %s!',fname); end

if exist('userLabels','var')
    if ischar(userLabels), userLabels={userLabels}; end  % must be cell of strings
    userLabels=unique(upper(userLabels));
end

%% Advance a single cell-by-cell record & compute the number or records in file

fprintf('\nTrying to read %s as BINARY file...',fname); 

bytes=contrRec(fp);

nbyt=ftell(fp); fseek(fp,0,'eof'); Nbyt=ftell(fp); nRec=Nbyt/nbyt; % number of records in file
if rem(Nbyt,nbyt),
   fprintf('failed\n');
   error('Budget file %s unknown binary format! Can''t help it, use a better compiler\n',fname);
end
fprintf('it works!\n');

%% Store record heading data, but do not yet the values,
% because we may perhaps only need a small selection

fprintf('Scanning %d headers\n',nRec);
B = repmat(struct('period', 0, 'tstp', 0,'label','','NLAY',0,'NROW',0,'NCOL',0),nRec, 1);
for i=1:nRec
    fseek(fp,(i-1)*nbyt,'bof');  % move to start of last layer in file
    [   B(i).period,...
        B(i).tstp,...
        B(i).label,...
        B(i).NLAY,...
        B(i).NROW,...
        B(i).NCOL]=contrRec(fp,bytes);
    if rem(i, 100)==0, fprintf('.'); end
    if rem(i,5000)==0, fprintf('%d records read\n',i); end
end
fprintf('finished, %d records read\n',i);
    
%% use indices into unique list of labels instead of labels themselves
% iLabInFile is index into the list of unique labels in the budget file

budlabels=unique({B.label}');
iLabInFile=NaN(nRec,1);
for i=1:nRec, iLabInFile(i)=strmatchi(B(i).label,budlabels); end

periodsInFile=[B.period]'; NPER = length(unique(periodsInFile));
tstpsInFile  =[B.tstp]';   NSTP = length(unique(tstpsInFile));
                           NLBL = length(budlabels);     % already unique
NCOL = B(end).NCOL;
NROW = B(end).NROW;
NLAY = B(end).NLAY;

pif=unique(periodsInFile);
tif=unique(tstpsInFile);

fprintf('File contains the following:\n');
fprintf('Stress Periods           : '); fprintf('% 3d',pif); fprintf('\n');
fprintf('Time steps               : '); fprintf('% 3d',tif); fprintf('\n');
fprintf('Number of records in file: %10d\n',nRec);
fprintf('Number of stress periods : %10d\n',length(pif));
fprintf('Number of time steps     : %10d\n',length(tif));
fprintf('Number of layers         : %10d\n',NLAY);
fprintf('Number of Rows           : %10d\n',NROW); 
fprintf('Number of columns        : %10d\n',NCOL);
fprintf('Number of unique labels  : %10d\n',NLBL);
for i=1:NLBL, fprintf('%s\n',budlabels{i}); end

if exist('userLabels','var') && isnumeric(userLabels) && userLabels(1)<=0,
    return;  % return the headers of all reords in file
end

%% select user specified selecton or all

fprintf('Please wait while I''m setting up selectors ...\n');


% translate userLabels into indices into list of unique labels
if exist('userLabels','var') && isempty(strmatch('',userLabels,'exact')) % verify if user specified non-empty labels
    % transl
    iUserLabels=zeros(size(userLabels));
    for iL=1:length(userLabels)
        iUserLabels(iL)=strmatchi(userLabels{iL},budlabels); % automatically generates an error if label illegal
    end
else % user want all
    userLabels=budlabels;              % all userLabels in the file in order
    iUserLabels=1:length(budlabels);   % their indices for storage
end


%% We create the same sort of selectors for the stress userPeriods
if  exist('userPeriods','var') && ~isempty(userPeriods)
    if ~isempty(userPeriods(userPeriods<1 | userPeriods>NPER))
        error('Requested userPeriods beyond range(1..%d) in budget file %s!',NPER,fname);
    end
    userPeriods=unique(userPeriods);
else
    userPeriods=unique(periodsInFile);
end

%% and also for the time steps
if exist('userTstps'  ,'var') && ~isempty(userTstps)
    if ~isempty(userTstps(userTstps<1 | userTstps>NSTP))
        error('Requested userPeriods beyond range(1..%d) in buget file %s!',NSTP,fname);
    end
    userTstps=unique(userTstps);
else
    userTstps=unique(tstpsInFile);
end

%% Selectors have values for all records of the inptut file. There
%  contents are userPeriod index, the userTstep index, the userLayers index
%  the combined userPeriod and userTstep index plus
%  the output record index and the output layer index
ULABC=1; UPCOL=2; UTCOL=3; UPTCOL=4;

Select=zeros(nRec,UPTCOL);  
for i=1:length(userLabels),  Select(iLabInFile    == iUserLabels(i),ULABC)=i; end
for i=1:length(userPeriods), Select(periodsInFile == userPeriods(i),UPCOL)=i; end
for i=1:length(userTstps),   Select(tstpsInFile   == userTstps(i)  ,UTCOL)=i; end

%% Rows and columns
if exist('userRows'   ,'var') && ~isempty(userRows)
    if ~isempty(userRows(userRows<1 | userRows>NROW))
        error('Requested userRows are beyond the range(1..%d) in budget file %s!',NROW,fname);
    end
else
    userRows=1:NROW;
end

if exist('userCols'   ,'var') && ~isempty(userCols)
    if ~isempty(userCols(userCols>NCOL | userCols<1))
        error('Requested columns are beyond the range(1..%d) in budget file %s!',NCOL,fname);
    end
else
    userCols=1:NCOL;
end

if exist('userLayers'   ,'var') && ~isempty(userLayers)
    if ~isempty(userLayers(userLayers<1 | userLayers>NLAY))
        error('Requested layers are beyond range(1..%d) in budget file %s!',NLAY,fname);
    end
else
    userLayers=(1:NLAY)'; % all layers
end

%% Get the unique period and tstep combinations for output records
%  records in inFile to deal with (LPT=Label Period TimeStep Layer
%  combinations)

%  = inFile records to deal with
LPTL=find(Select(:,ULABC)>0 & Select(:,UPCOL)>0 & Select(:,UTCOL)>0);

% unique period tstp combinations (PT combinations)
userPT = unique([periodsInFile(LPTL) tstpsInFile(LPTL)],'rows');

% output records in Select(:,UPTCOL)
for i=1:size(userPT,1),
    Select(periodsInFile ==userPT(i,1) & tstpsInFile==userPT(i,2),UPTCOL)=i;
end

nRecOut = length(unique(Select(LPTL,UPTCOL)));

if rem(nRecOut,1)~=0,
    error('mflab:readBud:nRecOutNotAnInteger',...
        ['nRecOut (number of output rectorcs = %g, is not an integer !\n',...
        'Check if MT3D or SEAWAT finished normally.\n',...
        'Check your LST file, to see if something went wrong during the simulation.'],nRecOut);
end


B=repmat(  struct(...
    'period' ,0,...
    'tstp'   ,0,...
    'label'  ,{[]},...
    'NLAY'   ,0,'NROW',0,'NCOL',0,...
    'term'   ,{[]},...
    'lays'   ,0,...
    'rows'   ,0,...
    'cols'   ,0),...
    nRecOut, 1);

%% Get the actual data values, here we go ...
fprintf('Reading requested data...\n');
for i=1:length(LPTL)
    iRecIn=LPTL(i);
    fseek(fp,nbyt*(iRecIn-1),'bof');
    iROut =Select(iRecIn,UPTCOL);   % output record Nr (sequential PT combination)
    iLabel=Select(iRecIn,ULABC);

    [B(iROut).period,...
        B(iROut).tstp,...
        label,...
        B(iROut).NLAY,...
        B(iROut).NROW,...
        B(iROut).NCOL,...
        values...
        ]=contrRec(fp,bytes);
    
    B(iROut).label{iLabel}=label;
    B(iROut).term{iLabel}=values(userRows,userCols,userLayers); % values is a 3D array in the Budget file
    B(iROut).lays=userLayers;
    B(iROut).rows=userRows;
    B(iROut).cols=userCols;

     if verbose>0
         fprintf('Label(period=%3d,tstp=%3d) =%s\n',B(iROut).period,B(iROut).tstp,B(iROut).label{iLabel});
     else
        if rem(i, 100)==0; fprintf('.');                      end
        if rem(i,5000)==0; fprintf('%d %d\n',iRecIn,iROut);   end
     end
end
fprintf('%6d records read.\n.',iRecIn);
fprintf('%6d records in output struct.\n',iROut);
fprintf('\n');

% That's all ... TO 090105

fclose(fp);
end

function [kper,kstp,label,nlay,nrow,ncol,values]=contrRec(fp,bytes)
% [kstp,kper,text,nlay,ncol,nrow,data]=contrRec(fid,bytes)
% --- reads a complete layer from a UNFORMATTED MODFLOW budget file
% --- if nargin==1, only the number of bytes at the beginning and enf
% of every record is returned. This is compiler dependent.
% once we know the number of bytes, a full record is returend.
% TO 070703 091214

if nargin==1, % just get the offset of this unformatted file
    fread(fp, 2,'int');  % kper, kstp
    label =char(fread(fp,16,'char')');
    
    % find offset from last null byte in label
    bytes=find(label==0,1,'last'); if isempty(bytes), bytes=0; end
    fread(fp,bytes,'int8'); % offset
    
    ncol=fread(fp, 1,'int');    % ncol nrow nlay
    nrow=fread(fp, 1,'int');    % ncol nrow nlay
    nlay=fread(fp, 1,'int');    % ncol nrow nlay
    
    fread(fp,bytes,'int8'); % offset
    fread(fp,bytes,'int8'); % offset
    
    n=ftell(fp); fread(fp,1,'float'); floatlen=ftell(fp)-n;
    fseek(fp,floatlen*(ncol*nrow*nlay-1),0);
    
    fread(fp,bytes,'int8'); % offset
    
    kper=bytes;
    return;
end

% offset is supposed to be give in bytes if narout>1
fread(fp,bytes,'int8');

kstp  =fread(fp, 1,'int');
kper  =fread(fp, 1,'int');
label =char(fread(fp,16,'char')');
label(label==' ')='';

ncol  =fread(fp, 1,'int');
nrow  =fread(fp, 1,'int');
nlay  =fread(fp, 1,'int');

fread(fp,bytes,'int8');
fread(fp,bytes,'int8');

if nargout<7, % don't fetch but skip
    n=ftell(fp); fread(fp,1,'float'); floatlen=ftell(fp)-n;
    fseek(fp,floatlen*(ncol*nrow*nlay-1),0);
else % fetch
    values=permute(reshape(fread(fp,ncol*nrow*nlay,'float'),[ncol,nrow,nlay]),[2,1,3]);  % is now nrow,ncol,nlay  (looking from above, Matlab style)    
end

fread(fp,bytes,'int8');

end
