function  [tsr,stat]=readTsrPoints(fname)
%READTSRPOINTS reads time series file produced by MODFPATH6
%
% Example:
%    [tsr,stat]=readTSR(basename,lpf);
%
% TO 130219

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

%0.
fprintf('# MATLAB %s %s\n',mfilename,datestr(now));

fid=fopen(fname,'r');
if fid<0, error('%s: Can''t open TSR file produced by MODPATH: <<%s>>',mfilename,fname); end

%% reading header

tsr.fname             = fname;
tsr.label             = fscanf(fid,'%s',1);
tsr.version           = fscanf(fid,'%d',1);
tsr.revision          = fscanf(fid,'%d',1);

%2
tsr.trackingDirection = fscanf(fid,'%d',1);
tsr.referenceTime     = fscanf(fid,'%f',1);

%% check if we are still on track:
         fscanf(fid,'%s',1);
header = fscanf(fid,'%s',1);
if ~strcmp(header,'HEADER')
    error(['%s: error reading tsr points file <<%s>>;\n',...
        'at this point we should have ''ENDHEADER'' instead of %s'],...
        mfilename,fname,header);
end

tsr.colHdr = {...
    'iTimePoint','iStp','time','id','iGrp',...
    'xG'  ,'yG'  ,'zG'  ,'iGrid',...
    'iLay','iRow','iCol',...
    'xL'  ,'yL'  ,'zL'};

%% getting the endpoints

fgets(fid); % move to next line

%% Determine the number of data lines in the file:
ip1  = ftell(fid);                  % remember position in file
fgets(fid);                         % read a full line including \b
ip2 = ftell(fid);                   % record start position in file
n     = ip2-ip1;                    % record length
fseek(fid,0,1); LFile = ftell(fid); % file size
nL    = (LFile - ip1)/n;            % number of data records in file
fseek(fid,ip1,-1);                  % reset file pointer to starting position
% check if the number of records is integer
if nL~=round(nL)
    error(['%s: modpath6 TSR file <<%s>> may be incomplete.\n',...
        'The number of point data records <<%d>> is not a whole number.'],...
        mfilename,fname,nL);
end

%% read the data backword to assue preallocation
tsr.P = NaN(nL,15);

for iLine = 1:nL
    tsr.P(iLine,:) = fscanf(fid,...
        ['%d %d %f %d %d',...
         ' %f %f %f %d',...
         ' %d %d %d',...
         ' %f %f %f'],[1,15]);
     fgets(fid);
end


%% Statistics if desired
tsr.nPoints     = size(tsr.P,1);
tsr.nTimePoints = numel(unique(tsr.P(:,1)));
tsr.nTimeSteps  = max(tsr.P(:,2));
tsr.tmin        = min(tsr.P(:,3));
tsr.tmax        = max(tsr.P(:,3));

fprintf('\n# Statistics for tsr file %s\n',tsr.fname);
fprintf('Total nr of points     = %d\n', tsr.nTimePoints);
fprintf('Total nr of time steps = %d\n', tsr.nTimeSteps);
fprintf('Min time               = %g\n',  tsr.tmin);
fprintf('Max time               = %g\n',  tsr.tmax);
fprintf('# ==============================\n');

if nargout == 2
    stat.header      = sprintf('%s version %d rivision %d\n',...
                        tsr.label,tsr.version,tsr.revision);
    stat.fname       = tsr.fname;
    stat.nTimePoints = tsr.nTimePoints;
    stat.nTimeSteps  = tsr.nTimeSteps;
    stat.tmin        = tsr.tmin;
    stat.tmax        = tsr.tmax;
end

fclose(fid);
