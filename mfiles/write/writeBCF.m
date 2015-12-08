function  writeBCF(basename,bcf)
%WRITEBCF writes input fiel for MODFLOW's basic flow package file (is BCF6)
%
% Example:
%    writeBCF(basename,bcf) --- 
%
% TO 0706030 081225 090713 091206

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if ~isfield(bcf,'FREE'), bcf.FREE = false; else bcf.FREE = logical(bcf.FREE); end

fid=fopen([basename,'.',bcf.ext],'wt');

%0. No header allowed in BCF file
 fprintf(    '%s\n',['# MODFLOW writeBCF ' datestr(now)]);

%1.
fprintf(fid,'%10d%10g%10d%10g%10d%10d     %s\n',...
    bcf.IBCFCB,bcf.HDRY,bcf.IWDFLG,bcf.WETFCT,bcf.IWETIT,bcf.IHDWET,...
    '     IBCFCB HDRY IWDFLG WETFCT IWETIT IHDWET');
%2
bcf.LAYCON(bcf.LAYCON<0)=0; bcf.LAYCON(bcf.LAYCON>3)=0;
bcf.LAYAVG(bcf.LAYAVG<0)=0; bcf.LAYAVG(bcf.LAYAVG>3)=0;

if bcf.FREE, fmt='  %1d%1d'; ibreak=20; else fmt='%1d%1d'; ibreak=40; end
 for iL=1:bcf.GRID.Nlay
     fprintf(fid,fmt,bcf.LAYAVG(iL),bcf.LAYCON(iL));
     if rem(iL,ibreak)==0 && iL<bcf.GRID.Nlay
         fprintf(fid,'\n');
     end
 end
fprintf(fid,'     LTYPE(NLAY) (40I2) LAYAVG|LAYCON\n');
%3 Anisotropy factor
warray(fid,bcf.TPRY,bcf.unit,'(10E15.5)','TPRY (hor anisotropy factor)',true,bcf.FREE);

%% computing VCONT required by the bcf package

PLANE=ones(bcf.GRID.Ny,bcf.GRID.Nx); % to change a layer par into a layer-wide cell par

k=0; iTRAN=0; iHY=0; % counters
for iL=1:bcf.GRID.Nlay
    k=k+1;
    %4 SF1
    if any(bcf.isTran)
        if numel(bcf.SF1)==1
            warray(fid,bcf.SF1        ,bcf.unit,'(10E15.5)',...
                sprintf('SF1(%d) LAYCON=%d',iL,bcf.LAYCON(iL)),true,bcf.FREE);
        else
            warray(fid,bcf.SF1(:,:,iL),bcf.unit,'(10E15.5)',...
                sprintf('SF1(%d) LAYCON=%d',iL,bcf.LAYCON(iL)),true,bcf.FREE);
        end
    end
    %5 TRAN where transmissiviy is fixed
    if any(bcf.LAYCON(iL)==[0 2])  % i.e. fixed transmissivity
        if numel(bcf.TRAN)==1
            warray(fid,bcf.TRAN        ,bcf.unit,'(10E15.5)',...
                sprintf('TRAN(%d) LAYCON=%d',iL,bcf.LAYCON(iL)),true,bcf.FREE);
        else
            iTRAN=iTRAN+1;
            warray(fid,bcf.TRAN(:,:,iTRAN),bcf.unit,'(10E15.5)',...
                sprintf('TRAN(%d) LAYCON=%d',iL,bcf.LAYCON(iL)),true,bcf.FREE);
        end
    %6 HY where layer is convertable
    else % Horizontal conductivity
        if ~isfield(bcf,'HY'),
            error(['%s: Need HY for kh in BCF package for convertible layers, i.e.: layer %d!\n', ...
                'Make sure HY in workspace agrees with LAYCON in worksheet LAY (BCF package)'],mfilename,iL);
        end
        if numel(bcf.HY)==1
            warray(fid,bcf.HY        ,bcf.unit,'(10E15.5)',...
                sprintf('HY(%d) LAYCON=%d',iL,bcf.LAYCON(iL)),true,bcf.FREE);
        else
            iHY=iHY+1;
            warray(fid,bcf.HY(:,:,iHY),bcf.unit,'(10E15.5)',...
                sprintf('HY(%d) LAYCON=%d',iL,bcf.LAYCON(iL)),true,bcf.FREE);
        end
    end
    %7 VCONT for all but the last layer
    if iL<bcf.GRID.Nlay  % write the vertical conductance between the two layers
        if numel(bcf.VCONT)==1
            warray(fid,bcf.VCONT        ,bcf.unit,'(10E15.5)',...
                sprintf('VCONT(%d) LAYCON=%d',iL,bcf.LAYCON(iL)),true,bcf.FREE);
        else
            warray(fid,bcf.VCONT(:,:,iL),bcf.unit,'(10E15.5)',...
                sprintf('VCONT(%d) LAYCON=%d',iL,bcf.LAYCON(iL)),true,bcf.FREE);
        end
    end
    %8 SF2 secondary storage coefficient = SY if LAYCON is 2 or 3
    if any(bcf.isTran) && any(bcf.LAYCON(iL)==[2 3])
        if ~isfield(bcf,'SF2'),
            error('Need SF2 for SY in BCF layer %d',iL);
        end 
        if numel(bcf.SF2)==1
            warray(fid,bcf.SF2        ,bcf.unit,'(10E15.5)',...
                sprintf('SF2(%d) LAYCON=%d',iL,bcf.LAYCON(iL)),true,bcf.FREE);
        else
            warray(fid,bcf.SF2(:,:,iL),bcf.unit,'(10E15.5)',...
                sprintf('SF2(%d) LAYCON=%d',iL,bcf.LAYCON(iL)),true,bcf.FREE);
        end
    end
    %9 WETDRY rewettabiliy only if LAYCON is 1 or 3
    if bcf.IWDFLG~=0 && any(bcf.LAYCON(iL)==[1 3])
        if iscell(bcf.WETDRY)
            if length(bcf.WETDRY{iL}(:))==1
                warray(fid,bcf.WETDRY{iL}*PLANE,bcf.unit,'(10E15.5)',sprintf(' WETDRY{%d} LAYCON=%d',iL,bcf.LAYCON(iL)),true,bcf.FREE);  %WETDRY
            else
                warray(fid,bcf.WETDRY{iL}      ,bcf.unit,'(10E15.5)',sprintf(' WETDRY{%d} LAYCON=%d',iL,bcf.LAYCON(iL)),true,bcf.FREE);  %WETDRY
            end
        else
            warray(fid,bcf.WETDRY(iL)*PLANE,bcf.unit,'(10E15.5)', sprintf(' Wetdry(%d) LAYCON=%d',iL,bcf.LAYCON(iL)),true,bcf.FREE);  %WETDRY
        end
    end

end
    
fclose(fid);
