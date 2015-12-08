function  huf=readHUF(fname,huf)
%READHUF reads Horizontal Unite Flow  package file
%
% Example:
%    huf=readHUF(basename,huf);
%
% TO 0706030 090713 090717

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

%0.
fprintf('# MATLAB readHUF %s\n',datestr(now));

fid=fopen(fname,'r');

fseek(fid,0,'eof'); eof=ftell(fid); fseek(fid,0,'bof');  % determine end of file

skipmodflowcomments(fid);

%1.
huf.IHUFCB=fscanf(fid,'%d',1);  % flag and budget file unit number
huf.HDRY  =fscanf(fid,'%f',1);  % head assigned to dry cells
huf.NHUF  =fscanf(fid,'%d',1);  % number of hydrologic units
huf.NPHUF =fscanf(fid,'%d',1);  % number of huf parameters
huf.IOHUF =fscanf(fid,'%d',1);  % >0 print calculated heads to unit IOHUF
fprintf(fgets(fid));

%2 Layer convertibility (confined / unconfined)
huf.LTHUF=mudread(fid,[huf.NLAY,1],'noctrlrec');

%3 Layer wet-dry switch 
huf.LAYWT=mudread(fid,[huf.NLAY,1],'noctrlrec');

%4  Defaults, don't touch, if at least one layer is wettable
if any(huf.LAYWT)
    huf.WETFCT=fscanf(fid,'%f',1);  % factor needed in calculation of head initially
    huf.IWETIT=fscanf(fid,'%d',1);  % iteration interval for attempting cell wetting
    huf.IHDWET=fscanf(fid,'%d',1);  % flag determining which equation to be used
    fgets(fid);
end

%5 WETDRY combination of layer threshold and flag indicating which
%neighboring cells can wet this cell
huf.WETDRY=zeros(huf.NROW,huf.NCOL,sum(huf.LAYWT>0));
k=0;
if any(huf.LAYWT)
    for i=1:huf.NLAY
        if huf.NLAY
            k=k+1;
            huf.WETDRY(:,:,k)=mudread(fid,[huf.NROW,huf.NCOL]);
        end
    end
end

%6
huf.TOP=zeros(huf.NROW,huf.NCOL,huf.NHUF);
for iHuf=1:huf.NHUF
    %6
    huf.HGUNAM{iHuf}=fscanf(fid,'%s',1); fgets(fid);
    %7
    huf.TOP(:,:,iHuf) =mudread(fid,[huf.NROW,huf.NCOL]);
    %8
    huf.THCK(:,:,iHuf)=mudread(fid,[huf.NROW,huf.NCOL]);
end

%9
for iHuf=1:huf.NHUF
    HGUNAM =fscanf(fid,'%s',1);
    HGUHANI=fscanf(fid,'%f',1);
    HGUVANI=fscanf(fid,'%f',1);
    fgets(fid);
    if strcmp(upper(HGUNAM),'ALL')
        for jHuf=1:huf.NHUF
            huf.HGUHANI(jHuf)=HGUHANI;
            huf.HGUVANI(jHif)=GHUVANI;
        end
        break;
    else
        found=0;
        for jHuf=1:huf.NHUF
            if strcmp(upper(HGUNAM),upper(huf.HGUNAM(jHuf)))
                huf.HGUHANI=HGUHANI;
                huf.HGUVANI=HGUVANI;
                found=1;
            end
        end
        if ~found
            error('HGUNAM=%s must be defined before reading HGUHANI and HGUVANI',HGUNAM);
        end
    end
end    

for iPar=1:huf.NPHUF
    huf.PARNAM{iPar}=fscanf(fid,'%s',1);
    huf.PARTYP{iPar}=fscanf(fid,'%s',1);
    huf.Parval(iPar)=fscanf(fid,'%f',1);
    huf.NCLU(iPar)  =fscanf(fid,'%d',1);
    fprintf(fgets(fid));
%9  parameter clusters
    for iClu=1:huf.NCLU(iPar)
        huf.Cluster(iPar).HGUNAM{iClu}=fscanf(fid,'%s',1);
        huf.Cluster(iPar).Mltarr{iClu}=fscanf(fid,'%s',1);
        huf.Cluster(iPar).Zonarr{iClu}=fscanf(fid,'%s',1);
        s=fgets(fid);
        huf.Cluster(iPar).IZ{iClu}=sscanf(s,'%d');
    end
end

%12, optional printcodes and printflags

huf.PRINTCODE=zeros(huf.NHUF,1);
huf.PRINTFLAGS=cell(huf.NHUF,1);

p=ftell(fid); % to properly read optional lines at the end of this file
while p<eof
    Prnt=fscanf(fid,'%s',1);
    HGUNAM=fscanf(fid,'%s',1);
    PRINTCODE=fscanf(fid,'%f',1);
    s=fgets(fid); PRINTFLAGS=sscanf(s,'%f',1);  % all flags on rest of line: HK HANI VK SS SY or ALL
    if ~strmp(upper(Prnt),'PRINT')
        error('Item 12 of HUF must start with the word ''PRINT'' not with %s\n',Print);
    end
    if strcmp(upper(HGUNAM),'ALL')
       for iHuf=1:huf.NHUF
           huf.PRINTCODE(iHuf)=PRINTCODE;
           huf.PRINTFLAGS{iHuf}=PRINTFLAGS;
       end
    else
        found=0;
        for jHuf=1:huf.NHUF
            if strcmp(upper(HGUNAM),upper(huf.HGUNAM{jHuf}))
                huf.PRINTCODE(jHuf)=PRINTCODE;
                huf.PRINTFLAGS{jHuf}=PRINTFLAGS;
                found=1;
                break;
            end
        end
        if ~found
            error('HGUNAM=%s must be defined before reading its printcodes',HGNUNAM);
        end
    end
    p=ftell(fid);
end

fclose(fid);
