function [Bgt,budlabels]=readBud(varargin)
%READBUD reads MODFLOW's budget file into a Matlab structure array.
%
% USAGE:
%   [B,budlabels]=readBud([+]fname [,userLabels [,userPeriods [,tsteps [,userLays [,userRows [,userCols]]]]]])
%
%   To add time, add 'time'time pair to arguments where time may be obtained from
%   heads: time = [H.totim]
%   To add pertim, add pair 'pertim',pertim to arguments where pertim may be
%   obtained from heads: pertim = [H.pertim]
%
%   ReadBud allows fast selection of any portion of the budget file, no
%   matter how big the file is.
%
%   budlabels are the list of unique labels in the entire budget file
%
% Example1: To get only the headers, i.e. meta data from the file, do this:
%    B=readBud(fname,-1);
%
% Examle2: to fetch all data from file do this
%    B=readBud(fname);
%
% Example3: to get more verbose screen output prepend fname with a "+" 
%    B=readBud(['+' fname]);
%
% Example4: 
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
% Example6: include time and pertim in the struct:
%    B = readBud(fname,'time',[H.totim],'pertim',[H.pertim]);
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
% TO 090104 091214
% TO 100602 make budlabels unique because labels may vary between stress userPeriods
% TO 100818 logic redone to clean up in accordance with readDat and readMT3D
% TO 100917 added labels to output
% TO 120205 major update to speed up selection, logic completely redone
% TO 130824 more flexible input using varargin
% SEE ALSO: readDat readMT3D readBud2
%
% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later
%

[time,       varargin] = getProp(varargin,{'t','totim'},[]);
[pertim,     varargin] = getProp(varargin,'pertim',[]);

[fname,      varargin] = getNext(varargin,'char',''); % always get filename first

[userLabels, varargin] = getNext(varargin,{'cell','char'},{});

[userPeriods,varargin] = getNext(varargin,'double',[]);
[userTstps,  varargin] = getNext(varargin,'double',[]);
[userLays,   varargin] = getNext(varargin,'double',[]);
[userRows,   varargin] = getNext(varargin,'double',[]);
[userCols,   ~       ] = getNext(varargin,'double',[]);


if     fname(1)=='+' || fname(1)=='1', fname=fname(2:end); verbose=1;
elseif fname(1)=='-' || fname(1)=='0', fname=fname(2:end); verbose=0;
else   verbose=0;
end

fp=fopen(fname);
if fp<1, error(['%s: Can''t find or open file %s!\n',...
        'REMEDY: In case you ran MT3DMS verify in the NAM sheet that you also\n',...
        '        ran one of the modflow codes MF2000, MF2005 ...'],mfilename,fname)
end

fseek(fp,0,1); n=ftell(fp); fseek(fp,0,-1);

if n==0
    error(['%s: file <<%s>> opened but has length <<%d>>.\n',...
        'Probably MODFLOW did not converge or finished early for some other\n',...
        'reason, like a crash due to reading bad or unreadible input.\n',...
        'When running MT3DMS, perhaps you did not switch on one of the MODFLOW programs\n',...
        'in the NAM worksheet. MT3DMS requires MODFLOW to run first.\n',...
        'REMEDY: check the list file for convergence and other reasons for an early halt\n',...
        '  with MT3DMS check one of the MODFLOW programs is on in the NAM sheet\n',...
        '  also verify that the package FTL and LMT6 are switched on in the NAM sheet.\n',...
        '  also check that any necesssary packages are on in the NAM sheet.'],...
        mfilename,fname,n);
end

doAllPeriods = isempty(userPeriods);
doAllTstps   = isempty(userTstps);
doAllLabels  = isempty(userLabels);
doAll        = doAllPeriods && doAllTstps && doAllLabels;

if ~doAllLabels && ~iscell(userLabels), userLabels={userLabels}; end

%% Advance a single cell-by-cell record & compute the number or records in file

fprintf('\nTrying to read %s as BINARY file...',fname); 

bytes=contrRec(fp);

nbyt=ftell(fp); fseek(fp,0,'eof'); Nbyt=ftell(fp); nRec=Nbyt/nbyt; % number of records in file
if rem(Nbyt,nbyt),
   fprintf('failed\n');
   error(['%s: Budget file <<%s>> has unknown binary format!\n',...
       'May be it''s incomplete or it is a COMPACT budget file.\n',...
       'Check the option OC/COMPACT in the MFLOW worksheet.\n',...
       'You may try to use readBud6 instead to read a compact budget file.\n'],mfilename,fname);
end
fprintf('it works!\n');

%% Store record heading data, but do not yet the values,
% because we may perhaps only need a small selection

% we will immediately check whether we need this record or not

fprintf('Scanning %d headers\n',nRec);

B = repmat(...              % allocate memory for headers of records in file
        struct('period', 0,...
        'tstp' , NaN,...
        'label','',...
        'NLAY' ,NaN,...
        'NROW' ,NaN,...
        'NCOL' ,NaN,...
        'mustSelect',0,...  % extract this record ?
        'iROut',NaN,...     % its output record number (compute)
        'iLbl' ,NaN),...    % its output label  number (compute)
    nRec, 1);

%% Get the header (all of them) and not the associated data

iRecOut=0; oldPeriod=-1; oldTStp=-1;

for i=1:nRec
    fseek(fp,(i-1)*nbyt,'bof');  % move to start of last layer in file
    [   B(i).period,...
        B(i).tstp,...
        B(i).label,...
        B(i).NLAY,...
        B(i).NROW,...
        B(i).NCOL]=contrRec(fp,bytes);

    %% Do we need this record, if so, set output record number
    % We have foreach(SP) --> timesteps --> foreach(timestep) --> labels
    
    newPeriod= B(i).period~=oldPeriod;

    newTStp  = B(i).tstp  ~=oldTStp;
    
    % start counting labels at each new time step
    if newPeriod || newTStp, iLbl=0; end

    B(i).mustSelect= doAll;
    
    if ~doAll
        B(i).mustSelect    = doAllPeriods || any(B(i).period==userPeriods);
        if B(i).mustSelect
           B(i).mustSelect = doAllTstps   || any(B(i).tstp==userTstps);
           if B(i).mustSelect
               B(i).mustSelect = doAllLabels  || strmatchi(B(i).label,userLabels,'noErr');
               if B(i).mustSelect
                   iLbl=iLbl+1;
               else
                   % B(i).mustSelect remains false
               end
           end
        end
    else
        iLbl=iLbl+1;
    end
    
    %if iLbl==0
    %    error('%s: Requested labels %s not found in budget file <<%s>>!',...
    %        mfilename,sprintfs('<<%s>>',userLabels),fname);
    %end
      
    if B(i).mustSelect,
        B(i).iLbl =iLbl;
        if newPeriod|| newTStp, iRecOut=iRecOut+1; end
        B(i).iROut=iRecOut;
        oldTStp   = B(i).tstp;
        oldPeriod = B(i).period;
    end
    
    if rem(i, 100)==0, fprintf('.'); end
    if rem(i,5000)==0, fprintf('%d records read\n',i); end
end
fprintf('finished, %d records scanned\n',i);
    
%% use indices into unique list of labels instead of labels themselves
% iLabInFile is index into the list of unique labels in the budget file

NCOL = B(end).NCOL;
NROW = B(end).NROW;
NLAY = B(end).NLAY;
NLBL = max([B.iLbl]);
NPER = max([B.period]);
NTSTP= max([B.tstp]);

fprintf('File contains the following:\n');
fprintf('Number of records in file: %10d\n',length(B));
fprintf('Number of stress periods : %10d\n',NPER);
fprintf('Number of time steps     : %10d\n',NTSTP);
fprintf('Number of layers         : %10d\n',NLAY);
fprintf('Number of Rows           : %10d\n',NROW); 
fprintf('Number of columns        : %10d\n',NCOL);
fprintf('Number of unique labels  : %10d\n',NLBL);
budlabels=unique({B.label});
for i=1:NLBL, fprintf('%s\n',budlabels{i}); end

%% Set the selected user Rows, Cols and Layers
if isempty(userRows), userRows=1:NROW; end
if isempty(userCols), userCols=1:NCOL; end
if isempty(userLays), userLays=1:NLAY; end

%% Getting the data
fprintf('\nReading the requested data ...\n');

nRecOut = max([B.iROut]);

if rem(nRecOut,1)~=0,
    error('mflab:readBud:nRecOutNotAnInteger',...
        ['nRecOut (number of output rectorcs = %g, is not an integer !\n',...
        'Check if MT3D or SEAWAT finished normally.\n',...
        'Check your LST file, to see if something went wrong during the simulation.'],nRecOut);
end


Bgt=repmat(...
    struct(...
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

fprintf('Please wait while I''m getting the requested data ...\n');

k=0;
for iRin = find([B.mustSelect]);
    
    fseek(fp,nbyt*(iRin-1),'bof');
    
    iROut= B(iRin).iROut;

    [Bgt(iROut).period,...
        Bgt(iROut).tstp,...
        label,...
        Bgt(iROut).NLAY,...
        Bgt(iROut).NROW,...
        Bgt(iROut).NCOL,...
        values...
        ]=contrRec(fp,bytes);
    
    iLbl = B(iRin).iLbl;
    
    Bgt(iROut).label{iLbl}=label;
    Bgt(iROut).term{ iLbl}=values(userRows,userCols,userLays); % values is a 3D array in the Budget file
    Bgt(iROut).lays=userLays;
    Bgt(iROut).rows=userRows;
    Bgt(iROut).cols=userCols;

    if verbose>0        
        fprintf('Label(period=%3d,tstp=%3d) =%s\n',...
            Bgt(iROut).period,Bgt(iROut).tstp,Bgt(iROut).label{iLabel});
    else        
        fprintf('.');
        k=k+1;
        if rem(k,50)==0; fprintf('%d\n',iROut);   end
    end
end

fprintf('%6d records in output struct.\n',length(Bgt));
fprintf('\n');

Iout = unique([B([B.mustSelect]).iROut]);
if any(diff(Iout)>1) % gaps ??
    msg ='mflab:readBud:missingOutrecs';
    warning('on',msg);
    cprintf('System','Missing stress period(s) in records read from %s: ',fname);
    cprintf('System',' %d',find(diff(Iout>1)));
    fprintf('\n');
    warning(msg,...
        ['%s: Records of some stress periods have not been read from <<%s>>.\n',...
        'This may ok, as, for example, not all records may have cell by cell flows that match your\n',...
        'requested label, as for instance wells do not have to be present at every stress period.\n',...
        'When using the output of readBud with this warning,\n',...
        'make sure you match the actual stress period and time step with\n,'...
        'that of other output variables such as heads and concentrations.'],...
        mfilename,fname);
    warning('off',msg);
end
% That's all ... TO 090105

fclose(fp);

% if time and or pertim were specified in the argument list, add them
if ~isempty(time) && numel(time)==numel(Bgt)
    for i=numel(Bgt):-1:1
        Bgt(i).time = time(i);
    end
end
if ~isempty(pertim) && numel(pertim) == numel(Bgt)
    for i=numel(Bgt):-1:1
        Bgt(i).pertim = pertim(i);
    end
end



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
    
    n=ftell(fp);
    fread(fp,1,'float');
    floatlen=ftell(fp)-n;
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
