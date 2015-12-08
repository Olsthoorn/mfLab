function  writeUPW(basename,upw)
%WRITEUPW writes input file for MODFLOW-NWT UPW pacakge (upwind flow)
%
% Example:
%    writeUPW(basename,uwp) --- write linear flow package file
%
% TO 130327

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

writeLPF(basename,upw);
