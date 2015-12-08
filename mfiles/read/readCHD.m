function chd=readCHD(fname,pth,chd)
%READCHD reads MODFLOW's CHD boundary package file
%
% Example:
%    chd=readCHD(fname,pth,chd) 
%
% TO 070711 090714

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

%0
fprintf('# MATLAB readCHD %s\n',datestr(now));
fid=fopen([pth fname],'r');
skipmodflowcomments(fid);

% 1 [PARAMETER NPCHD MXL]
s=fgets(fid); C=textscan(s,'%s %d %d',1);
if strcmp(upper(C{1}{1}),'PARAMETER'),
    fprintf(s);
    chd.NPCHD=C{2};
    chd.MXL  =C{3};
    s=fgets(fid);
else
    chd.NPCHD=0;
    chd.MXL  =0;
end

%2
fprintf(s);
C=textscan(s,'%s %s',1);  % May be there is an error in the manual as output budget units misses in input instructions
chd.MXACTC=sscanf(C{1}{1},'%d',1); % max number of active chders during any stress period
chd.ICHDCB=sscanf(C{2}{1},'%d',1); % unit or writing wel buget values 

%% Read auxiliary parameters if present
NAux=size(C{1},1)-1;
for iAux=1:NAux
    chd.Aux{iAux}=C{2}{iAux+1};
end

%% format for cell values
fmt=' %f %f %f %f %f %f %f %f %f %f'; fmt=fmt(1:3*(5+NAux));

%% 3 Parameters
for iPar=1:chd.NPCHD
    s=fgets(fid); fprintf(s);
    C=textscan(s,'%s %s %f %d',1); 
    chd.PARNAM(iPar)=C{1};
    chd.PARTYP(iPar)=C{2};
    chd.PARVAL(iPar)=C{3};
    chd.NLST(iPar)  =C{4};

    %4
    for inlst=1:chd.NLST(iPar)
        s=fgets(fid); fprintf(s);
        C=textscan(s,fmt,1,'CollectOutput',1);
        chd.par(iPar).values(inlst,1:5+NAux)=C{1};
    end
end

%% Stress periods

for iPer=1:chd.NPER
    %5
    s=fgets(fid); fprintf(s);
    C=textscan(s,'%d %d',1);
    chd.ITMP(iPer)=C{1};  % Number of non-parameter well data to be read (-1 is reuse)
    chd.NP(iPer)  =C{2};  % Number of parameters in use in current period
    
    %6
    if chd.ITMP(iPer)>0
        fprintf('reading chd\n');
        tic;
        chd.cel(iPer).values=textscan(fid,fmt,chd.ITMP);
        toc;
        fgets(fid);
    end
    
    %7
    if chd.NP(iPer)>0
        for i=1:chd.NP(iPer)
            s=fgets(fid); fprintf(s);
            C=textscan(s,'%s',1);
            chd.Pname(iPer).name(i)=C{1};
        end
    end
end

fclose(fid);
