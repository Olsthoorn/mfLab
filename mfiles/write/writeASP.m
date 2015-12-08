function writeASP(basename,asp)
%WRITEASP writes input file for Doherty's MODFLOW-ASP program
%
% Example:
%    writeASP(basename,asp)
%
% See also: MODFLOW-NWT
%
% TO 100830

% Copyright 2010-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fid=fopen([basename,'.',asp.ext],'wt');

%0.
fprintf('# MATLAB  writeASP %s\n',datestr(now));
fprintf('The only functionality is a better rewetting than mf2k is capable of.\n');
fprintf('All other more pest related options are unused till now.\n');

asp.IPESTINT=0; asp.INTERP=0;  % to only suse ASP for better rewetting functionality

%1. value used for inactive cells
fprintf(fid,'%10d%10d     IPESTINT INTERP\n',asp.IPESTINT,asp.INTERP);

%2. value used for inactive cells
fprintf(fid,'%10d%15g%10d%15g     NOSTOP HDRYBOT LIMOP MINTHICK\n',...
    asp.NOSTOP,asp.HDRYBOT,asp.LIMOP,asp.MINTHICK);

fclose(fid);
