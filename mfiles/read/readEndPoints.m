function  [endp,stat]=readEndPoints(fname)
%READENDPOINTS read simulation end points produced by MODPATH6
%
% Example:
%    [endp,stat]=readEndp(basename,lpf);
%
% TO 130219

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

% The column headers of the data array:
endp.colHdr = {...
    'id','iGrp','iStatus',...
    'tI','tF',...
    'iGridI','iLayI','iRowI','iColI','iFaceI','iZoneI',...
    'xLI','yLI','zLI','xGI','yGI','zGI',...
    'iGridF','iLayF','iRowF','iColF','iFaceF','iZoneF',...
    'xLF','yLF','zLF','xGF','yGF','zGF','label'};

nFld = numel(endp.colHdr)-1;  % -1 because of label of different type

%0.
fprintf('# MATLAB %s %s\n',mfilename,datestr(now));
  
%%
fid=fopen(fname,'r');
if fid<0, error('%s: Can''t open Endp file produced by MODPATH: <<%s>>',mfilename,fname); end

fseek(fid,0,1);
if ftell(fid)==0
    error('file %s is empty',fname);
    return;
end
fseek(fid,0,-1);

%% reading header

endp.fname             = fname;
endp.label             = fscanf(fid,'%s',1);
endp.version           = fscanf(fid,'%d',1);
endp.revision          = fscanf(fid,'%d',1);

%2
endp.trackingDirection = fscanf(fid,'%d',1);
endp.totalCount        = fscanf(fid,'%d',1);
endp.releaseCount      = fscanf(fid,'%d',1);
endp.maximumID         = fscanf(fid,'%d',1);
endp.referenceTime     = fscanf(fid,'%f',1);

%3
endp.pending           = fscanf(fid,'%d',1);
endp.active            = fscanf(fid,'%d',1);
endp.normallyTerminated= fscanf(fid,'%d',1);
endp.zoneTerminated    = fscanf(fid,'%d',1);
endp.unreleased        = fscanf(fid,'%d',1);
endp.stranded          = fscanf(fid,'%d',1);
%4
endp.groupCount        = fscanf(fid,'%d',1);
endp.groupNm{endp.groupCount} = 'groupName';
%5
for iGrp=1:endp.groupCount
    endp.groupNm{iGrp}=   fscanf(fid,'%s',1);
    fgets(fid); % proceed to start of new line
end

%% check if we are still on track:
         fscanf(fid,'%s',1);
header = fscanf(fid,'%s',1);
if ~strcmp('HEADER',header)
    error(['%s: error reading end points file <<%s>>;\n',...
        'at this point we should have ''ENDHEADER'' instead of %s'],...
        mfilename,fname,header);
end

fgets(fid);

%% Determine the number of data lines in the file:
% ip1  = ftell(fid);                  % remember position in file
% fgets(fid);                         % read a full line including \b
% ip2 = ftell(fid);                   % record start position in file
% n     = ip2-ip1;                    % record length
% fseek(fid,0,1); LFile = ftell(fid); % file size
% nL    = (LFile - ip1)/n;            % number of data records in file
% fseek(fid,ip1,-1);                  % reset file pointer to starting position
% % check if the number of records is integer
% if nL~=round(nL)
% end

%% read the data backword to assue preallocation
endp.P = NaN(endp.totalCount,nFld);
endp.Lbl{endp.totalCount,1} = 'dummy';

try
    for iLine = 1:endp.totalCount
        endp.P(iLine,:) = fscanf(fid,...
            ['%d %d %d %f %f',...
             ' %d %d %d %d %d %d',...
             ' %f %f %f %f %f %f',...
             ' %d %d %d %d %d %d',...
             ' %f %f %f %f %f %f'],[1,nFld]);
         endp.Lbl{iLine} = fscanf(fid,'%s',1);
         fgets(fid);
    end
catch ME
    fprintf(1,ME.message); fprintf(1,'\n');
    msgId = 'mfLab:mfPath:readEndPoints_incomplete_file';
    warning(msgId,'on');
    warning(msgId,...
        ['%s: modpath6 Endp file <<%s>> may be incomplete.\n',...
        'The number of point data records <<%d>> is not a whole number.'],...
        mfilename,fname,iLine,endp.totalCount);
    warning(msgId,'off');
end

%% Statistics if desired
if nargout==2
    stat.nPoints     = size(endp.P,1);
    stat.tmin        = min(endp.P(:,3));
    stat.tmax        = max(endp.P(:,3));

    stat.fname             = endp.fname;
    stat.label             = endp.label;
    stat.version           = endp.version;
    stat.revision          = endp.revision;

    %2
    stat.trackingDirection = endp.trackingDirection;
    stat.totalCount        = endp.totalCount;
    stat.releaseCount      = endp.releaseCount;
    stat.maximumID         = endp.maximumID;
    stat.referenceTime     = endp.referenceTime;

    %3
    stat.pending           = endp.pending;
    stat.active            = endp.active;
    stat.normallyTerminated= endp.normallyTerminated;
    stat.zoneTerminated    = endp.zoneTerminated;
    stat.unreleased        = endp.unreleased;
    stat.stranded          = endp.stranded;
    %4
    stat.groupCount        = endp.groupCount;
    %5
    stat.groupNm           = endp.groupNm;

    fprintf('\n# Statistics for endp file %s\n',endp.fname);
    fprintf('%s version %d revision\n',endp.label,endp.version,endp.revision);

    fldNames = fieldnames(stat);

    for ifld = 1:numel(fldNames)
        aField = fldNames{ifld};

        fprintf('%-20s = ',aField);
        switch class(stat.(aField))
            case 'char'
                fprintf('%s\n',stat.(aField));
            case 'double'
                fprintf('%g\n',stat.(aField));
            case 'cell'
                for j=1:numel(stat.(aField))
                    fprintf(' "%s"',stat.(aField){j});
                end
                fprintf('\n');
            otherwise
                error('%s: unexpected class <<%s>> for field <<%s>> in stat',...
                    mfilename,class(stat.(aField)),aField);
        end
    end
    fprintf('# ==============================\n');
end

fclose(fid);
