function Obs=readOBS(fname)
%READOBS reads formatted conc obsrvation point outut 9(MT3D, Seawat)
%
% Example:
%    Obs=readObs(fname);
%
% TO 091012

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fid=fopen(fname); if fid<1, error('READOBS: Can''t find or open file %s!',fname); end


fprintf(fgets(fid));

Obs.Loc=sscanf(fgets(fid),'%d %d %d',[3,Inf])'; % trasnpose
Obs.NOBS=size(Obs.Loc,1);

fmt=['%d %f', repmat(' %g',1,Obs.NOBS)];
data=fscanf(fid,fmt,[2+Obs.NOBS,Inf])';  % transpose
Obs.trstp=data(:,1);
Obs.t    =data(:,2);
Obs.values=data(:,3:end);

fclose(fid);
