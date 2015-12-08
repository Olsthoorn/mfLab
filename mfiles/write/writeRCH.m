function writeRCH(basename,rch)
%WRITERCH writes input file for MODFLOW's recharge (RCH) package
%
% Example:
%    writeRCH(basename,rch);
%
% TO 091007

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fid=fopen([basename,'.RCH'],'wt');

%0.
fprintf(fid,'# MATLAB writeRCH, %s\n',datestr(now));
fprintf(    '# MATLAB writeRCH, %s\n',datestr(now));

%1 no parameters use in recharge
fprintf(fid,'PARAMETER     0\n');

%2. NRCHOP IRCHCB
% NRCHOP
%    1=recharge only in toplayer of model
%    2=recharge into specific layer specified in IRCH
%    3=recharge into highest active layer
% IRCHCB  unit number of budgetfile to save budgetdata
fprintf(fid,'%10i%10i  NRCHOP IRCHCB\n',[rch.NRCHOP,rch.IRCHCB]);

%3 for each stress period
for iP=1:rch.NPER
    %5 INRECH INIRCH
    %  INRECH
    %   >=0 layer variable of recharge values is used
    % in mfLab if INRECH==0 then RECH must be defined in workspace
    %   <0  use recharge from previous stress period
    %  INIRCH
    %   >0  # of parameters used to define RECH
    %   <0  recharge parameters of previous stress period are used
    %    0  no parameters (not allowed, omit??)
    fprintf(fid,'%10d%10d  INRECH(%d)  INIRCH(%d)\n',...
        rch.INRECH(iP),rch.INIRCH(iP),iP,iP);
    
    %6
    if rch.INRECH(iP)>=0
        if iscell(rch.RECH)
            warray(fid,rch.RECH{iP},rch.unit,'(10E12.3)',sprintf(' RECH{%d}',iP),true,rch.FREE);
        else
            warray(fid,rch.RECH(:,:,iP),rch.unit,'(10E12.3)',sprintf(' RECH{%d}',iP),true,rch.FREE);
        end
    end
    
    %7   parameterized items not used
    %8
    if rch.NRCHOP==2 && rch.INIRCH(iP)>=0
        if iscell(rch.IRCH)
            warray(fid,rch.IRCH{iP},rch.unit,'(20I4)',sprintf(' IRCH{%d}',iP),true,rch.FREE);
        else
            warray(fid,rch.IRCH(:,:,iP),rch.unit,'(20I4)',sprintf(' IRCH{%d}',iP),true,rch.FREE);
        end
    end
end
fclose(fid);
