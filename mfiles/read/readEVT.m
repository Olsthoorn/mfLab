function evt=readEVT(fname,evt)
%READEVT reads modflow's EVT package input file
%
% Example:
%    evt=readEVT(fname,evt);
%
% TO 090714

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

%0.
fprintf('# MATLAB readEVT %s\n',datestr(now));

fid=fopen(fname,'r');
skipmodflowcomments(fid);

% 1 [PARAMETER NPEVT]
s=fgets(fid); C=textscan(s,'%s %d',1);
if strcmp(C{1}{1},'PARAMETER'),
    fprintf(s);
    evt.NPEVT=C{2};
    s=fgets(fid);
else
    evt.NPEVT=0;
end

%2. NEVTOP IEVTCB
% NEVTOP
%    1=recharge only in toplayer of model
%    2=recharge into specific layer specified in IEVT
%    3=recharge into highest active layer
% IEVTCB  unit number of budgetfile to save budgetdata
fprintf(s);
C=textscan(s,'%s %s',1);
evt.NEVTOP=sscanf(C{1}{1},'%d',1); % max number of active wells during any stress period
evt.IEVTCB=sscanf(C{2}{1},'%d',1); % unit or writing wel buget values 

%3 Parameters
for iPar=1:evt.NPEVT
    s=fgets(fid); fprintf(s);
    C=textscan(s,'%s %s %f %d',1); 
    evt.PARNAM(iPar)=C{1};
    evt.PARTYP(iPar)=C{2};
    evt.PARVAL(iPar)=C{3};
    evt.NCLU(iPar)  =C{4};
    evt.Cluster(iPar)={[]};
   
    %4
    for iClu=1:evt.NCLU(iPar);
        evt.Cluster{iPar}.Mltarr{iClu}=fscanf(fid,'%s',1);
        Zonarr=fscanf(fid,'%s',1);
        evt.Cluster{iPar}.Zonarr{iClu}=Zonarr;
        s=fgets(fid);
        if strcmp(Zonarr,'ALL')
            fprintf(s);
        else
            evt.Cluster{iPar}.iZ=sscanf(s,'%d');
        end
    end
end

% FOR EACH STRESS PERIOD
for iP=1:evt.NPER
    s=fgets(fid); fprintf(s);
    C=textscan(s,'%d %d %d %d',1);
    %5 INRECH INIEVT
    evt.INSURF(iP)=C{1};
    evt.INEVTR(iP)=C{2};
    evt.INEXDP(iP)=C{3};
    if evt.NEVTOP==2
        evt.INIEVT(iP)=C{4};
    end

    %6
    if evt.INSURF(iP)>=0
        evt.SURF{iP}=mudread(fid,[evt.NROW,evt.NCOL]);
    end
    
    %7   parameters Pname [IEVTPF]
    if evt.NPEVT==0 && evt.INEVTR(iP)>=0
        evt.EVTR{iP}=mudread(fid,[evt.NROW,evt.NCOL]);
    end
       
    %8
    if evt.NPEVT>0 && evt.INEVTR>0
        for inr=1:evt.INEVTR(iP)
            s=fgets(fid); fprintf(s);
            C=textscan(s,'%s %d',1);
            
            % Pname
            evt.Param{iPar}.Pname{inr}=C{1};
            
            % [IEVTPF]
            if ~isempty(C{2})
                evt.Param{iPar}.IEVTPF{inr}=C{2};
            end
        end
    end

    %9
    if evt.INEXDP>=0
        evt.EXDP{iP}=mudread(fid,[evt.NROW,evt.NCOL]);
    end
        
    %10
    if evt.NEVTOP==2 && evt.INIEVT(iP)>=0
        evt.IEVT{iP}=mudread(fid,[evt.NROW,evt.NCOL]);
    end
end
fclose(fid);
