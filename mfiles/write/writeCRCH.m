function writeCRCH(basename,crch)
%WRITECRCH writes input file for conduit flow recharge package
%
% Example:
%    writeCRCH(basename,crch)
%
% TO 090708 090713 090715 110513

% Copyright 2009-2011 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fprintf('# MATLAB writeCRCH %s\n',datestr(now));

fid=fopen([basename,'.',crch.ext],'wt');

% Item 0-1 
fprintf(fid,'# item 0: Conduit recharge file\n');
fprintf(fid,'%10d\n',crch.IFLAG_CRCH);

% Item 2
fprintf(fid,'%10d %10f\n',[crch.NODE_NUMBERS crch.P_CRCH]');

fclose(fid);
