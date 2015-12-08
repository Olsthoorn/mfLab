function  lpf=readLPF(fname,lpf)
%READLPF reads layer property flow package file
%
% Example:
% lpf=readLPF(basename,lpf);
%
% TO 0706030 090713

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

%0.
fprintf('# MATLAB readLPF %s\n',datestr(now));

fid=fopen(fname,'r');

skipmodflowcomments(fid);

%1.
lpf.ILPFCB=fscanf(fid,'%d',1);
lpf.HDRY  =fscanf(fid,'%f',1);
lpf.NPLPF =fscanf(fid,'%d',1);
fprintf(fgets(fid));

%2 Layer convertibility (confined / unconfined)
lpf.LAYTYP=mudread(fid,[lpf.NLAY,1],'noctrlrec');

%3 Interblock conductance computation method
lpf.LAYAVG=mudread(fid,[lpf.NLAY,1],'noctrlrec');

%4 Horizontal layer anisotropy
lpf.CHANI=mudread(fid,[lpf.NLAY,1],'noctrlrec');

%5 vertical anisotropy (we always use 1 in the spreadsheeet
lpf.LAYVKA=mudread(fid,[lpf.NLAY,1],'noctrlrec');

%6 Layer wet-dry switch 
lpf.LAYWET=mudread(fid,[lpf.NLAY,1],'noctrlrec');

%7  Defaults, don't touch, if at least one layer is wettable
if any(lpf.LAYWET)
    lpf.WETFCT=fscanf(fid,'%f',1);
    lpf.IWETIT=fscanf(fid,'%d',1);
    lpf.IHDWET=fscanf(fid,'%d',1);
    fgets(fid);
end

%8  we don't use parameters
if lpf.NPLPF>0
%    N=lpf.NPLPF;
%    lpf.PARNAM=cell(N,1);
%    lpf.PARTYP=cell(N,1);
%    lpf.parval=NaN(N,1);
%    lpf.NCLU  =NaN(N,1);
%    lpf.Cluster=cell(N,1);
    
    for iPar=1:lpf.NPLPF
        lpf.PARNAM{iPar}=fscanf(fid,'%s',1);
        lpf.PARTYP{iPar}=fscanf(fid,'%s',1);
        lpf.parval(iPar)=fscanf(fid,'%f',1);
        lpf.NCLU(iPar)  =fscanf(fid,'%d',1);
        fprintf(fgets(fid));
    %9  parameter clusters
        for iClu=1:lpf.NCLU(iPar)
            lpf.Cluster{iPar}.layer(iClu) =fscanf(fid,'%d',1);
            lpf.Cluster{iPar}.Mltarr{iClu}=fscanf(fid,'%s',1);
            lpf.Cluster{iPar}.Zonarr{iClu}=fscanf(fid,'%s',1);
            s=fgets(fid);
            lpf.Cluster{iPar}.IZ{iClu}    =fscanf(s,'%d');
        end
    end
end

%% Allocate memory
lpf.KH=zeros(lpf.NROW,lpf.NCOL,lpf.NLAY);
lpf.KV=zeros(lpf.NROW,lpf.NCOL,lpf.NLAY);
if any(lpf.LAYWET) && any(lpf.LAYCON)
    lpf.WETDRY=zeros(lpf.NCOL,lpf.NROW,lpf.NLAY);
end
if lpf.isTran
    lpf.SS=zeros(lpf.NROW,lpf.NCOL,lpf.NLAY);
    if any(lpf.LAYTYP)
        lpf.SY=zeros(lpf.NROW,lpf.NCOL,lpf.NLAY);
    end
end

%% Get the layer data
for iL=1:lpf.NLAY
    %10.
    if lpf.NPLPF>0 && strmatchi('HK',lpf.PARTYP,'ErrOpt')>0
        s=fgets(fid); lpf.HK_Printcode(iL)=sscanf(s,'%d',1);
    else
        lpf.KH(:,:,iL)=mudread(fid,[lpf.NROW,lpf.NCOL]);  % hiorzontal conductivity
    end

    %11 hor anisotropy through CHANI not through item 11
     if lpf.CHANI(iL)<=0
        if lpf.NPLPF>0 && strmatchi('HANI',lpf.PARTYP,'ErrOpt')>0
            s=fgets(fid); lpf.HANI_Printcode(iL)=sscanf(s,'%d',1);
        else
            lpf.HANI(:,:,iL)=mudread(fid,[lpf.NCOL,lpf.NROW]);
        end
     end
    
    %12 vertical conductivity
    if lpf.NPLPF>0 && (strmatchi('VK',lpf.PARTYP,'EorrOpt')>0 || strmatchi('VANI',lpf.PARTYP,'ErrOpt')>0)
        s=fgets(fid); lpf.VKA_Printcode(iL)=sscanf(s,'%d',1);
    else
        lpf.KV(:,:,iL)=mudread(fid,[lpf.NROW,lpf.NCOL]);  % vertical conductivity    
    end
    
    if lpf.isTran
        
        %13 Ss
        if lpf.NPLPF>0 && strmatchi('SS',lpf.PARTYP,'EorrOpt')>0
            s=fgets(fid); lpf.SS_Printcode(iL)=sscanf(s,'%d',1);
        else
            lpf.SS(:,:,iL)=mudread(fid,[lpf.NROW,lpf.NCOL]);  % Ss
        end
        
        %14 Sy
        if lpf.LAYTYP(iL)
            if lpf.NPLPF>0 && strmatchi('SY',lpf.PARTYP,'EorrOpt')>0
               s=fgets(fid); lpf.SY_Printcode(iL)=sscanf(s,'%d',1);
            else    
                lpf.SY(:,:,iL)=mudread(fid,[lpf.NROW,lpf.NCOL]);  %Sy
            end
        end
    end
    
    %15 skip is c for resitance layer below model layers
    if lpf.LAYCBD(iL)
        if lpf.NPLPF>0 && strmatchi('VKCB',lpf.PARTYP,'EorrOpt')>0
           s=fgets(fid); lpf.SY_Printcode(iL)=sscanf(s,'%d',1);
        else    
            lpf.VKCB(:,:,iL)=mudread(fid,[lpf.NROW,lpf.NCOL]);  %VKCB
        end
    end

    %16
    if lpf.LAYWET(iL) && lpf.LAYTYP(iL)
        lpf.WETDRY(:,:,iL)=mudread(fid,[lpf.NROW,lpf.NCOL]);  %WETDRY
    end
end    

fclose(fid);
