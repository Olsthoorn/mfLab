function bcn=readBCN(fname,bcn,type)
%READBCN reads MODFLOW's stress files (Boundary Condition Files)
%
% USAGE:
%    readBCN(fname,pth,bcn,type);
%
% type is one of WEL,DRN,RIV,GHB,CHD
%
% TO  090831 120323

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

%0
fprintf('# MATLAB readBCN for %s; %s\n',type,datestr(now));

fid=fopen(fname,'r');
skipmodflowcomments(fid);

%% 1 [PARAMETER NPRIV MXL]
bcn.NPar =0;
bcn.MXL  =0;

p=ftell(fid);
s=fgets(fid); C=textscan(s,'%s %d %d',1);
if strcmpi(C{1}{1},'PARAMETER'),
    %fprintf(s);
    bcn.NPar=C{2};
    bcn.MXL =C{3};
else
    fseek(fid,p,-1); % set pointer back
end

%% 2
bcn.MXACTD=fscanf(fid,'%10d',1);
s=fgets(fid);
%bcn.type=sscanf(s,'%d',1);  % optional, needs not be present
bcn.type=type;

%% 3
for iPar=1:bcn.NPar
    bcn.PARNAM(iPar)=fscanf(fid,'%s',1);
    bcn.PARTYP(iPar)=fscanf(fid,'%s',1);
    bcn.Parval(iPar)=fscanf(fid,'%d',1);
    bcn.NLST  (iPar)=fscanf(fid,'%d',1);
    fgets(fid);
    
    % look ahead to see how many values are on the next line
    p=ftell(fp); s=fgets(fp); n=length(scanf(s,'%f',[1,Inf])); fseek(fid,p,-1);
    fmt=repmat('%10f',n);
    bcn.parvals=textscan(fid,fmt,bcn.NLST(iPar));

        %4
    for inlst=1:bcn.NLST(iPar)
        error('Parameter reading not yet implemented in readBCN\n');
        %s=fgets(fid); fprintf(s);
        %C=textscan(s,fmt,1,'CollectOutput',1);
        %bcn.par(iPar).values(inlst,1:6+NAux)=C{1};
    end
 end

%% The non parameter values for this boundary condition
for iPer=1:bcn.NPER   %length(bcn.ITMP)
    s=fgets(fid); % fprintf(s);
    C=textscan(s,'%d %d',1);
    bcn.ITMP(iPer)=C{1};  % Number of non-parameter well data to be read (-1 is reuse)
    bcn.NP(iPer)  =C{2};  % Number of parameters in use in current period
    
    %6
    if bcn.ITMP(iPer)>0
        % look ahead to see how many values are on the next line
        p=ftell(fid); s=fgets(fid); n=length(sscanf(s,'%f',[1,Inf])); fseek(fid,p,-1);
        fmt=repmat('%10f',[1,n]);
        
        %then read bcn.ITMP(iPer) lines with this number of values
        % this assumes that there are no comments
        % else we have to adopt the procedure by adding info on the number
        % of fields for each specific boundary condition, which is a bit
        % risky because some packages are changed somewhat over time, like
        % the desnity options for CHD included with seawat version 4.
        % As long as we don't stumble over comments, we can use this simple
        % routine to read in all boundary conditions except multi-node
        % wells.
        bcn.cel(iPer)=textscan(fid,fmt,bcn.ITMP(iPer),'CollectOutput',1);
        %toc;
        fgets(fid);
    end

    %7
    if bcn.NP(iPer)>0
        for i=1:bcn.NP(iPer)
            bcn.Pname(iPer).name(i)=textscan(s,'%s %*[^\n]',1);
        end
    end
end

fclose(fid);
