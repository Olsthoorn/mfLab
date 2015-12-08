function swi=readSWI(fname,swi)
%READSWI reads Salt Water Intrusion input file
%
% Example:
%    swi=readSWI(fname,swi)
%
% TO 090831

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fid=fopen(fname,'r');

%% 0
fprintf('# MATLAB readSWI %s\n',datestr(now));

%% 1

swi.NPLN  =fscanf(fid,'%10d',1);
swi.ISTRAT=fscanf(fid,'%10d',1);
swi.ISWIZT=fscanf(fid,'%10d',1);
swi.NPRN  =fscanf(fid,'%10d',1);
fgets(fid);
    
%% 2 
swi.TOESLOPE=fscanf(fid,'%10f',1);
swi.TIPSLOPE=fscanf(fid,'%10f',1);
swi.ZETAMIN =fscanf(fid,'%10f',1);
swi.DELZETA =fscanf(fid,'%10f',1);
fgets(fid);

%% 3
% Values for dimensionless density (NPLN+1 if ISTRAT=1, NPLN if ISTRT=0
if swi.ISTRAT==1, n=swi.NPLN+1; else n=swi.NPLN+2; end
swi.NU=mudread(fid,[1,n]);

%% 4
% SWI needs all interfaces within all layers. So we intersect each interaces with the top and
% bottom of each layer on the fly and write out the results immediately according to the
% order required by SWI.

swi.zeta{swi.NPLN}=NaN(swi.NROW,swi.NCOL,swi.NLAY); % allocate
for ipln=1:swi.NPLN  % for each surface
    for iLay=1:swi.NLAY % for each layer
        swi.zeta{ipln}(:,:,iLay)=mudread(fid,[swi.NROW,swi.NCOL]);
    end
end

%% 5
% Types the value of the effective porosity by layer
swi.SSW=NaN(swi.NROW,swi.NCOL,swi.NLAY);
for iLay=1:swi.NLAY
    swi.SSW(:,:,iLay)=mudread(fid,[swi.NROW,swi.NCOL]);
end

%% 6
% Types the kind of source at each point by layer
swi.ISOURCE=NaN(swi.NROW,swi.NCOL,swi.NLAY);
for iLay=1:swi.NLAY
    swi.ISOURCE(:,:,iLay)=mudread(fid,[swi.NROW,swi.NCOL]);
end

%% The .swi file is ready to be used
fclose(fid);
