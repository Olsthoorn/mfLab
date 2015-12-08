function rch=readRCH(fname,rch)
%READRCH reads MODFLOW's recharge package input file (RCH)
%
% Example:
%    rch=readRCH(fname,rch);
%
% TO 090101

%0.
fprintf('# MATLAB readRCH %s\n',datestr(now));

fid=fopen(fname,'r');
skipmodflowcomments(fid);

% 1 [PARAMETER NPRCH]
s=fgets(fid); C=textscan(s,'%s %d',1);
if strcmpi(C{1}{1},'PARAMETER'),
    %fprintf(s);
    rch.NPRCH=C{2};
    s=fgets(fid);
else
    rch.NPRCH=0;
end

%2. NRCHOP IRCHCB
% NRCHOP
%    1=recharge only in toplayer of model
%    2=recharge into specific layer specified in IRCH
%    3=recharge into highest active layer
% IRCHCB  unit number of budgetfile to save budgetdata
%fprintf(s);
C=textscan(s,'%s %s',1);
rch.NRCHOP=sscanf(C{1}{1},'%d',1); % max number of active wells during any stress period
rch.IRCHCB=sscanf(C{2}{1},'%d',1); % unit or writing wel buget values 

%3 Parameters
for iPar=1:rch.NPRCH
    s=fgets(fid); %fprintf(s);
    C=textscan(s,'%s %s %f %d',1); 
    rch.PARNAM(iPar)=C{1};
    rch.PARTYP(iPar)=C{2};
    rch.PARVAL(iPar)=C{3};
    rch.NCLU(iPar)  =C{4};
    rch.Cluster(iPar)={[]};
   
    %4
    for iClu=1:rch.NCLU(iPar);
        rch.Cluster{iPar}.Mltarr{iClu}=fscanf(fid,'%s',1);
        Zonarr=fscanf(fid,'%s',1);
        rch.Cluster{iPar}.Zonarr{iClu}=Zonarr;
        s=fgets(fid);
        if strcmpi(Zonarr,'ALL')
            %fprintf(s);
        else
            rch.Cluster{iPar}.iZ=sscanf(s,'%d');
        end
    end
end

% FOR EACH STRESS PERIOD
for iP=1:rch.NPER
    s=fgets(fid); % fprintf(s);
    C=textscan(s,'%d %d',1);
    %5 INRECH INIRCH
    %  INRECH
    %   >=0 layer variable of recharge values is used
    %   <0  use recharge from previous stress period
    rch.INRECH(iP)=C{1};

    %5  INIRCH
    %   >0  # of parameters used to define RECH
    %   <0  recharge parameters of previous stress period are used
    %    0  no parameters (not allowed, omit??)
    if rch.NRCHOP==2,
        rch.INIRCH(iP)=C{2};
    end
    
    %6
    if rch.NPRCH==0 && rch.INRECH(iP)>=0
        rch.RECH{iP}=mudread(fid,[rch.NROW,rch.NCOL]);
    end
    
    %7   parameters Pname [IRCHPF]
    if rch.NPRCH>0 && rch.INRECH(iP)>0
        for inr=1:rch.INRECH(iP)
            s=fgets(fid); % fprintf(s);
            C=textscan(s,'%s %d',1);
            
            % Pname
            rch.Param{iPar}.Pname{inr}=C{1};
            
            % [IRCHPF]
            if ~isempty(C{2})
                rch.Param{iPar}.IRCHPF{inr}=C{2};
            end
        end
    end
        
    %8
    if rch.NRCHOP==2 && rch.INIRCH(iP)>=0
        rch.INRICH{iP}=mudread(fid,[rch.NROW,rch.NCOL]);
    end
end
fclose(fid);
