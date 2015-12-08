function moc=readMOCMAIN(fname,pth)
%READMOCMAIN reads mocDense main file
%
% Example:
%    moc=readMOCMAIN(fname,pth,moc);
%
% TO 090814, old BAS necessary for reading MOCDENSE


% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

%0.
fprintf('# MATLAB readMOCMAIN (MocDense main file) %s\n',datestr(now));

if pth(end)==filesep,
    fid=fopen([pth fname],'r');
else
    fid=fopen([pth filesep fname],'r');
end
%1 two heading lines
fprintf('%s',fgets(fid));
%2
fprintf('%s',fgets(fid));

%3.  NLAY NROW NCOL NPER ITMUNI
fscanf(fid,'%10d',1); moc.NLAY =fscanf(fid,'%10d',1);
fscanf(fid,'%10d',1); moc.NROW =fscanf(fid,'%10d',1);
fscanf(fid,'%10d',1); moc.NCOL =fscanf(fid,'%10d',1);
fgets(fid);

%% ======== SKIP 5 lines ====================================
fgets(fid); %Don't know wat this is, just skip the whole line
fgets(fid); %Don't know wat this is, just skip the whole line
fgets(fid); %Don't know wat this is, just skip the whole line
fgets(fid); %Don't know wat this is, just skip the whole line
fgets(fid); %Don't know wat this is, just skip the whole line
%% \======== SKIP 5 lines ===================================

%6.  if <0 constant head if > 0 compute head if 0 inactive (or lake in lake package)
moc.CHLOR=NaN(moc.NROW,moc.NCOL,moc.NLAY);
for ilay=1:moc.NLAY
    moc.CHLOR(:,:,ilay)=mudread(fid,[moc.NROW,moc.NCOL]);
end

%% ======== SKIP 5 lines ====================================
fgets(fid); %Don't know wat this is, just skip the whole line
fgets(fid); %Don't know wat this is, just skip the whole line
fgets(fid); %Don't know wat this is, just skip the whole line
fgets(fid); %Don't know wat this is, just skip the whole line
fgets(fid); %Don't know wat this is, just skip the whole line
%% \======== SKIP 5 lines ===================================

moc.SFLAG=NaN(moc.NROW,moc.NCOL,moc.NLAY);
for ilay=1:moc.NLAY
    moc.SFLAG(:,:,ilay)=mudread(fid,[moc.NROW,moc.NCOL]);
end

moc.DISPL =mudread(fid,[1,moc.NLAY]);
moc.DISPHT=mudread(fid,[1,moc.NLAY]);
moc.DISPVT=mudread(fid,[1,moc.NLAY]);
moc.RETARD=mudread(fid,[1,moc.NLAY]);

for ilay=1:moc.NLAY
    moc.DZ(  :,:,ilay)=mudread(fid,[moc.NROW,moc.NCOL]);
    moc.PEFF(:,:,ilay)=mudread(fid,[moc.NROW,moc.NCOL]);
end

fclose(fid);
