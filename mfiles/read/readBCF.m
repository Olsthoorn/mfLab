function  bcf=readBCF(fname,bcf)
%READBCF reads MODFLOW's basic flow package file, old version
%
% Example:
%    readBCF(basename,bcf) --- 
%
% TO 090814


% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fid=fopen(fname,'r');

fprintf('# MATLAB readBCF (old version) %s\n',datestr(now));

%1
bcf.ISS=fscanf(fid,'%10d',1);  % steady state flag if 0 then transient 1 then steady state
bcf.IBCFCB =fscanf(fid,'%d',1);
fgets(fid);

bcf.isTran=~bcf.ISS;

%2
bcf.LAYCON=fscanf(fid,'%2d',bcf.NLAY); fgets(fid);

%3 Anisotropy factor
bcf.TPRY=mudread(fid,[bcf.NLAY,1]);

%4
bcf.DELC=mudread(fid,[1,bcf.NCOL]);

%5
bcf.DELR=mudread(fid,[bcf.NROW,1]);

%% VCONT required by the bcf package

for i=1:bcf.NLAY
    %6 SF1
    if bcf.isTran
        if  bcf.LAYCON(i)==1
            bcf.SY=mudread(fid,[bcf.NROW,bcf.NCOL]);    % Primary storage   layer %d,  SF1=SY';
        else % LAYCON=0, 2 or 3 then SF1 is elastic aquifer storage coefficient
            bcf.SS=mudread(fid,[bcf.NROW,bcf.NCOL]);    % Primary storage   layer %d,  SF1=SS*DZ'
        end
    end
    %7 TRAN
    if bcf.LAYCON(i)==0  || bcf.LAYCON(i)==2
        bcf.TRAN(:,:,i)=mudread(fid,[bcf.NROW,bcf.NCOL]); % Transmissivity    layer %d, (along rows)
    else
    %8 KH (horizontal K (HY in manual of BCF, HK in BCF6)
        bcf.KH( :,:,i)  =mudread(fid,[bcf.NROW,bcf.NCOL]);   % Hydr. Conduct. KH layer %d, (along rows)
    %9 KH
        bcf.BOT(:,:,i)  =mudread(fid,[bcf.NROW,bcf.NCOL]);   % Bottom of layer %d, (along rows)
    end
    %10 VCONT
    if i<bcf.NLAY  % write the vertical conductance between the two layers
        bcf.VCONT(:,:,i)=mudread(fid,[bcf.NROW,bcf.NCOL]);  % Leakance VCONT
    end
    %8 SF2
    if bcf.isTran && (bcf.LAYCON(i)==2 || bcf.LAYCON(i)==3)
        bcf.SY(:,:,i)=mudread(fid,[bcf.NROW,bcf.NCOL]); % Secondary Storage layer %d, SF2=SY',1));
    end
    if bcf.LAYCON(i)==2 || bcf.LAYCON(i)==3
        bcf.TOP(:,:,i)=mudread(fid,[bcf.NROW,bcf.NCOL]); % Secondary Storage layer %d, SF2=SY',1));
    end   
end
    
fclose(fid);
