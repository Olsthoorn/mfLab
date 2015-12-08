function dis=readDIS(fname)
%READDIS reads MODFLOW's discretization package input file
%
% Example:
%    readDIS(basename,dis);
%
% TO 070630 090713

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fprintf('# MATLAB readDIS %s\n',datestr(now));

fid=fopen(fname,'r');

skipmodflowcomments(fid)

%1.
s=fgets(fid); % fprintf(s);
C=textscan(s,'%f %f %f %f %f %f',1);
dis.NLAY  =C{1};
dis.NROW  =C{2};
dis.NCOL  =C{3};
dis.NPER  =C{4};
dis.ITMUNI=C{5};
dis.LENUNI=C{6};

%2  Resistance at bottom of layer?? -- Quasi 3D flag for each layer
dis.LAYCBD=mudread(fid,[dis.NLAY,1],'norec');
dis.LAYCBD(end)=0; % never a LAYCBD layer below the lowest aquifer

%3.
dis.DELR=mudread(fid,[1,dis.NCOL]);

%4.
dis.DELC=mudread(fid,[dis.NROW,1]);


%% 5, reading top and bottom of all layers
dis.Z=NaN(dis.NROW,dis.NCOL,dis.NLAY); % allocate

% top of model
dis.Z(:,:,1)=mudread(fid,[dis.NROW,dis.NCOL]); %TOP of model');

%6 bottom of all layers
k=1;
for i=1:dis.NLAY
    k=k+1;
    dis.Z(:,:,k)=mudread(fid,[dis.NROW,dis.NCOL]);  % bottom of layer
    if dis.LAYCBD(i) && i<dis.NLAY
        k=k+1;
        dis.Z(:,:,k)=mudread(fid,[dis.NROW,dis.NCOL]); % bottom of LAYBCD;
    end
end

%% FOR EACH STRESS PERIOD
%7. 
dis.PERLEN=NaN(dis.NPER,1);
dis.NSTP  =NaN(dis.NPER,1);
dis.TSMULT=NaN(dis.NPER,1);
dis.TrSS  =cell(dis.NPER,1);
dis.isTran=NaN(dis.NPER,1);

for i=1:dis.NPER
    s=fgets(fid); % fprintf(s);
    C=textscan(s,'%f %d %f %s',1);
    dis.PERLEN(i)=C{1};
    dis.NSTP(i)  =C{2};
    dis.TSMULT(i)=C{3};
    dis.TrSS{i}  =C{4};
    dis.isTran(i)=strcmp('TR',dis.TrSS{i});
end

fclose(fid);
