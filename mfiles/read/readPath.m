function  [pth,stat]=readPath(fname)
%READPATH reads path line file produced by MODPATH6
%
% Example:
%    [pth,stat]=readPATHPOINTS(basename,lpf)
%
% TO 130219

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

%0.
fprintf('# MATLAB %s %s\n',mfilename,datestr(now));

fid=fopen(fname,'r');
if fid<0, error('%s: Can''t open path line file produced by MODPATH: <<%s>>',mfilename,fname); end

%% reading header

pth.fname             = fname;
pth.label             = fscanf(fid,'%s',1);
pth.version           = fscanf(fid,'%d',1);
pth.revision          = fscanf(fid,'%d',1);

%2
pth.trackingDirection = fscanf(fid,'%d',1);
pth.referenceTime     = fscanf(fid,'%f',1);

%% check if we are still on track:
         fscanf(fid,'%s',1);
header = fscanf(fid,'%s',1);
if ~strcmp(header,'HEADER')
    error(['%s: error reading pth points file <<%s>>;\n',...
        'at this point we should have ''ENDHEADER'' instead of %s'],...
        mfilename,fname,header);
end

pth.colHdr = {...
    'id','iGrp','iTimePoint','iStp','time',...
    'xG'  ,'yG'  ,'zG'  ,...
    'iLay','iRow','iCol','iGrid',...
    'xL'  ,'yL'  ,'zL','iLine'};

Nfld = numel(pth.colHdr);

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
pth.P = NaN(nL,Nfld);

fprintf('\n%s: Reading <<%d>> records from path line file <<%s>> ...\n',mfilename,nL,fname);

for iLine = 1:nL
    pth.P(iLine,:) = fscanf(fid,...
        ['%d %d %d %d %f',...
         ' %f %f %f',...
         ' %d %d %d %d',...
         ' %f %f %f %d'],[1,Nfld]);
     fgets(fid);
     if rem(iLine, 1000)==0, fprintf('.'); end
     if rem(iLine,50000)==0, fprintf('%d\n',iLine); end
end

fprintf('...\n <<%d>> lines read (thanks for waiting)\n',iLine);

%% Statistics if desired
pth.nPoints     = size(pth.P,1);
pth.nTimePoints = numel(unique(pth.P(:,3)));
pth.nTimeSteps  = max(pth.P(:,4));
pth.tmin        = min(pth.P(:,5));
pth.tmax        = max(pth.P(:,5));

fprintf('\n# Statistics for path file %s\n',pth.fname);
fprintf('Total nr of points     = %d\n', pth.nTimePoints);
fprintf('Total nr of time steps = %d\n', pth.nTimeSteps);
fprintf('Min time               = %g\n',  pth.tmin);
fprintf('Max time               = %g\n',  pth.tmax);
fprintf('# ==============================\n');

if nargout == 2
    stat.header      = sprintf('%s version %d rivision %d\n',...
                        pth.label,pth.version,pth.revision);
    stat.fname       = pth.fname;
    stat.nTimePoints = pth.nTimePoints;
    stat.nTimeSteps  = pth.nTimeSteps;
    stat.tmin        = pth.tmin;
    stat.tmax        = pth.tmax;
end

fclose(fid);
