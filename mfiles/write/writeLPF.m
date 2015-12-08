function  writeLPF(basename,lpf)
%WRITELPF writes input file for MODFLOW's LPF package (Layer Property Flow)
%
% Example:
%    writeLPF(basename,lpf) --- write linear flow package file
%
% TO 0706030

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fid=fopen([basename,'.',lpf.ext],'wt');

ctrRecDesired = true;
%0.
fprintf(fid,'# MATLAB  write%s %s\n',lpf.name,datestr(now));
fprintf(    '# MATLAB  write%s %s\n',lpf.name,datestr(now));

%1.
fprintf(fid,'%10d%10.3g%10d',lpf.I___CB,lpf.HDRY,lpf.NP___);

if lpf.upw, fprintf(fid,'%10d',lpf.IPHDRY); end

% MF2005 options
if ~lpf.upw
    if isfield(lpf,'STORAGECOEFFICIENT'), fprintf(fid,'  STORAGECOEFFICIENT'); end
    if isfield(lpf,'CONSTANTCV'),         fprintf(fid,'  CONSTANTCV'); end
    if isfield(lpf,'THICKSTRT'),          fprintf(fid,'  THICKSTRT'); end
    if isfield(lpf,'NOCVCORRECTION'),     fprintf(fid,'  NOCVCORRECTION'); end
end

fprintf(fid,'     I%sCB HDRY NP%s IPHDRY\n',lpf.name,lpf.name);

%2 Layer convertibility (confined / unconfined)
warray(fid,lpf.LAYCON(:)',lpf.unit,'(25I3)','LAYCON',~ctrRecDesired);

%3 Interblock conductance computation method
warray(fid,lpf.LAYAVG(:)',lpf.unit,'(25I3)','LAYAVG',~ctrRecDesired);

%4 Horizontal layer anisotropy
warray(fid,lpf.CHANI(:)',lpf.unit,'(10G10.3)','CHANI',~ctrRecDesired);

%5 vertical anisotropy: preferably use 0 indicating we provide VK directly
warray(fid,lpf.LAYVKA(:)',lpf.unit,'(25I3)','LAYVKA',~ctrRecDesired);

%6 Layer wet-dry switch 
warray(fid,lpf.LAYWET(:)',lpf.unit,'(25I3)','LAYWET',~ctrRecDesired);

%7  Layer wetting options
if any(lpf.LAYWET)
    fprintf(fid,'%10.3g%10d%10d     WETFCT IWETIT IHDWET\n',lpf.WETFCT,lpf.IWETIT,lpf.IHDWET);
end

%8  we don't use parameters
%9  we don't use parameters

PLANE=ones(lpf.GRID.Ny,lpf.GRID.Nx); % to change a layer par into a layer-wide cell par

transient=any(lpf.isTran);

iLAYCBD=0;

for iL=1:lpf.GRID.Nlay,
    %10.
    if numel(lpf.HK)==1
        warray(fid,lpf.HK        ,lpf.unit,'(10E15.5)',sprintf(' HK{%d}',iL),ctrRecDesired,lpf.FREE);  % horizontal conductivity
    else
        warray(fid,lpf.HK(:,:,iL),lpf.unit,'(10E15.5)',sprintf(' HK{%d}',iL),ctrRecDesired,lpf.FREE);  % horizontal conductivity
    end
    %11 hor anisotropy through CHANI not through item 11
    if lpf.CHANI(iL)<=0
        if numel(lpf.HANI)==1
            warray(fid,lpf.HANI        ,lpf.unit,'(10E15.5)',sprintf(' Hor anisotropy{%d}',iL),ctrRecDesired,lpf.FREE);  % horizontal anisotropy
        else
            warray(fid,lpf.HANI(:,:,iL),lpf.unit,'(10E15.5)',sprintf(' Hor anisotropy{%d}',iL),ctrRecDesired,lpf.FREE);  % horizontal anisotropy
        end
    end           
    %12 vertical conductivity
    if lpf.LAYVKA(iL)==0,
        str=sprintf(' VK{%d}',iL);
    else
        str=sprintf(' Vert Anisotropy{%d}',iL);
    end
    if numel(lpf.VKA)==1
        warray(fid,lpf.VKA        ,lpf.unit,'(10E15.5)',str,ctrRecDesired,lpf.FREE);  % vertical conductivity
    else
        warray(fid,lpf.VKA(:,:,iL),lpf.unit,'(10E15.5)',str,ctrRecDesired,lpf.FREE);  % vertical conductivity
    end
    if transient
        %13 Ss
        if ~isfield(lpf,'SS')
           error('SS was not specified in mf_adapt for transient simulation');
        end
        if numel(lpf.SS)==1
            warray(fid,lpf.SS        ,lpf.unit,'(10E15.5)',sprintf(' SS{%d}',iL),ctrRecDesired,lpf.FREE);  % Ss
        else
            warray(fid,lpf.SS(:,:,iL),lpf.unit,'(10E15.5)',sprintf(' SS{%d}',iL),ctrRecDesired,lpf.FREE);  % Ss
        end
        if lpf.LAYCON(iL)>0
            %14 Sy
            if ~isfield(lpf,'SY')
               error('SY was not specified in mf_adapt for transient simulation for converible layer %d',iL);
            end
            if numel(lpf.SY)==1
                warray(fid,lpf.SY        ,lpf.unit,'(10E15.5)',sprintf(' SY{%d}',iL),ctrRecDesired,lpf.FREE);  %Sy
            else
                warray(fid,lpf.SY(:,:,iL),lpf.unit,'(10E15.5)',sprintf(' SY{%d}',iL),ctrRecDesired,lpf.FREE);  %Sy
            end
        end
    end
    
    %15 Resistance layer below model layers
    if lpf.GRID.LAYCBD(iL) && iL~=lpf.GRID.Nlay
        
        iLAYCBD=iLAYCBD+1;
        
        if numel(lpf.VKCB)==1
            warray(fid,lpf.VKCB             ,lpf.unit,'(10E15.5)',sprintf(' VKCB{%d}',iL),ctrRecDesired,lpf.FREE);  %VKCB
        else
            warray(fid,lpf.VKCB(:,:,iLAYCBD),lpf.unit,'(10E15.5)',sprintf(' VKCB{%d}',iL),ctrRecDesired,lpf.FREE);  %VKCB
        end
    end
    
    %16
    if lpf.LAYWET(iL)&& lpf.LAYCON(iL);
        if iscell(lpf.WETDRY)
            if length(lpf.WETDRY{iL}(:))==1
                warray(fid,lpf.WETDRY{iL}*PLANE,lpf.unit,'(10E15.5)',sprintf(' WETDRY{%d}',iL),ctrRecDesired,lpf.FREE);  %WETDRY
            else
                warray(fid,lpf.WETDRY{iL}      ,lpf.unit,'(10E15.5)',sprintf(' WETDRY{%d}',iL),ctrRecDesired,lpf.FREE);  %WETDRY
            end
        else
            warray(fid,    lpf.WETDRY(iL)*PLANE,lpf.unit,'(10E15.5)',sprintf(' WETDRY{%d}',iL),ctrRecDesired,lpf.FREE);  %WETDRY
        end
    end
end    

fclose(fid);
