function  bcf=readBC6(fname,bcf)
%READBC6 reads MODLFOW's basic flow package file
%
% Example:
%    readBCF(basename,bc6);
%
% TO 0706030 081225 090713

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fid=fopen(fname,'r');

fprintf('# MATLAB readBCF6 %s\n',datestr(now));

%0
skipmodflowcomments(fid);

%1.
    bcf.IBCFCB =fscanf(fid,'%d',1);
    bcf.HDRY   =fscanf(fid,'%f',1);
    bcf.IWDFLG =fscanf(fid,'%d',1);
    bcf.WETFCT =fscanf(fid,'%f',1);
    bcf.IWETIT =fscanf(fid,'%d',1);
    bcf.IHDWET =fscanf(fid,'%d',1);
    fgets(fid);

%2
bcf.LTYPE=fscanf(fid,'%2d',bcf.NLAY); fgets(fid);
bcf.LAYCON=rem(bcf.LTYPE,10);
bcf.LAYAVG=(bcf.LTYPE-bcf.LAYCON)/10;

%3 Anisotropy factor
bcf.TPRY=mudread(fid,[bcf.NLAY,1]);

%% VCONT required by the bcf package

for i=1:bcf.NLAY
    %4 SF1
    if bcf.isTran
        if  bcf.LAYCON(i)==1
            bcf.SY=mudread(fid,[bcf.NROW,bcf.NCOL]);    % Primary storage   layer %d,  SF1=SY';
        else % LYACON=0, 2 or 3 then SF1 is elastic aquifer storage coefficient
            bcf.SS=mudread(fid,[bcf.NROW,bcf.NCOL]);    % Primary storage   layer %d,  SF1=SS*DZ'
        end
    end
    %5 TRAN
    if bcf.LAYCON(i)==0  || bcf.LAYCON(i)==2
        bcf.TRAN(:,:,i)=mudread(fid,[bcf.NROW,bcf.NCOL]); % Transmissivity    layer %d, (along rows)
    %6 KH
    else
        bcf.KH(:,:,i)  =mudread(fid,[bcf.NROW,bcf.NCOL]);   % Hydr. Conduct. KH layer %d, (along rows)
    end
    %7 VCONT
    if i<bcf.NLAY  % write the vertical conductance between the two layers
        bcf.VCONT(:,:,i)=mudread(fid,[bcf.NROW,bcf.NCOL]);  % Leakance VCONT
    end
    %8 SF2
    if bcf.isTran && (bcf.LAYCON(i)==2 || bcf.LAYCON(i)==3)
        bcf.SY(:,:,i)=mudread(fid,[bcf.NROW,bcf.NCOL]); % Secondary Storage layer %d, SF2=SY',1));
    end
    %9 WETDRY
    if bcf.IWDFLG~=0 && (bcf.LAYCON(i)==1 || bcf.LAYCON(i)==3)
        bcf.WETDRY(:,:,i)=mudread(fid,[bcf.NROW,bcf.NCOL]); % Secondary Storage layer %d, SF2=SY',1));
    end
    
end
    
fclose(fid);
