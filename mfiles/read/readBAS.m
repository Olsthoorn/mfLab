function bas=readBAS(fname,bas)
%READBAS read MODFLOW's basic bas6 package file
%
% Example:
%    bas=readBAS6(fname,pth,bas)
%
% TO 090814, old BAS necessary for reading MOCDENSE

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

%0.
fprintf('# MATLAB readBAS (Old version) %s\n',datestr(now));

fid=fopen(fname,'r');

%1 two heading lines
fprintf('%s',fgets(fid));
%2
fprintf('%s',fgets(fid));

%3.  NLAY NROW NCOL NPER ITMUNI
bas.NLAY =fscanf(fid,'%10d',1);
bas.NROW =fscanf(fid,'%10d',1);
bas.NCOL =fscanf(fid,'%10d',1);
bas.NPER =fscanf(fid,'%10d',1);
bas.ITMUNI=fscanf(fid,'%10d',1);
fgets(fid);

%4 (BCF WEL DRN RIV EVT XXX GHB RCH SIP XXX SOR OC)
s=fgets(fid);
bas.word=sscanf(s,'%s',1);
if ~strcmpi(bas.word,'FREE')
    bas.ITUNIT=sscanf(s,'%3d',[1,12]);
end

%5
bas.IAPART=fscanf(fid,'%10d',1);
bas.ISTRT =fscanf(fid,'%10d',1);
fgets(fid);

%6.  if <0 constant head if > 0 compute head if 0 inactive (or lake in lake package)
bas.IBOUND=NaN(bas.NROW,bas.NCOL,bas.NLAY);
for ilay=1:bas.NLAY
    bas.IBOUND(:,:,ilay)=mudread(fid,[bas.NROW,bas.NCOL]);
end

%3. value used for inactive cells
bas.HNOFLO=fscanf(fid,'%f',1); fprintf(fgets(fid));

%4. initial heads
bas.STRTHD=NaN(bas.NROW,bas.NCOL,bas.NLAY);
for ilay=1:bas.NLAY
    bas.STRTHD(:,:,ilay)=mudread(fid,[bas.NROW,bas.NCOL]);
end

fclose(fid);
