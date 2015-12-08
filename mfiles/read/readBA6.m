function bas=readBA6(fname,bas)
%READBA6 reads MODFLOW's bas6 package file
%
% Example:
%     bas=readBAS6(fname,bas);
%
% TO 070630, 090713 090714

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

%0.
fprintf('# MATLAB readBAS6 %s\n',datestr(now));

fid=fopen(fname,'r');
skipmodflowcomments(fid);

%1.  optional words XSECTION CHTOCH FREE PRINTTIME SHOWPROGRESS
bas.OPTIONS=fgets(fid);

%2.  if <0 constant head if > 0 compute head if 0 inactive (or lake in lake package)
bas.IBOUND=NaN(bas.NROW,bas.NCOL,bas.NLAY);
for ilay=1:bas.NLAY
    bas.IBOUND(:,:,ilay)=rarray(fid,[bas.NROW,bas.NCOL]);
end

%3. value used for inactive cells
bas.HNOFLO=fscanf(fid,'%f',1);
fgets(fid);

%4. initial heads
bas.STRTHD=NaN(bas.NROW,bas.NCOL,bas.NLAY);
for ilay=1:bas.NLAY
    bas.STRTHD(:,:,ilay)=rarray(fid,[bas.NROW,bas.NCOL]);
end

fclose(fid);
