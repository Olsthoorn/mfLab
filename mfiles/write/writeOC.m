function writeOC(basename,oc)
%WRITEOC writes input file for MODFLOW's output control file (OC)
%
% Example:
%    writeOC(fname,oc) --- write output control file
%
% TO 070630 081226 100911
 
% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fid=fopen([basename, '.',oc.ext],'wt');

%0
fprintf(fid,'# MATLAB writeOC %s\n',datestr(now));
fprintf(    '# MATLAB writeOC %s\n',datestr(now));

%1 IHEDFM IDDNFM IHEDUN IDDNUN (For each simulation)
fprintf(fid,'%10d%10d%10d%10d     IHEDFM IDDNFM IHEDUN IDDNUN \n',...
    oc.IHEDFM,oc.IDDNFM,oc.IHEDUN,oc.IDDNUN);
    %  IHEDFM head     print format code 0..12
    %  IDDNFM drawdown print format code 0..12
    %  IHEDUN unit number of file to save heads
    %  IDDNUN unit number of flle to save drawdowns


% % first oc.INCODE=0 to make sure that at least one following record is
% read specifying what to print
% and make sure all values are <1 so that always all layers are treated the
% same, mfLab does not implement layer-specific output
oc.INCODE(oc.INCODE>0)=0;
oc.INCODE(1)=0;

%2
for iP=1:oc.NPER  % for each stress period

    k=0; % in-stress-period time-step counter; reset at beginning of each stress period
    
    for iS=1:oc.NSTP(iP)  % for all time steps within this stress period

        k=k+1; % counting of total number of time steps
        
        ihddfl=oc.IHDDFL(iP) & (rem(k,oc.IHDDFL(iP))==0 | k==oc.NSTP(iP));
        ibudfl=oc.IBUDFL(iP) & (rem(k,oc.IBUDFL(iP))==0 | k==oc.NSTP(iP));
        icbcfl=oc.ICBCFL(iP) & (rem(k,oc.ICBCFL(iP))==0 | k==oc.NSTP(iP));
        
        %2 INCODE IHDDFL IBUDSV ICBCFL
        fprintf(fid,'%10d%10d%10d%10d     INCODE(%d) IHDDFL(%d) IBUDFL(%d) ICBCFL(%d)\n',...
            oc.INCODE(iP),ihddfl,ibudfl,icbcfl,k,k,k,k);
        % IHHDFL
        %   0= don't print heads or drawdown whatever item 3 says
        %   <>0 use item 3 codes to decide on printing
        %   and use abs(IHHDFL) as print frequency
        % IBUDFL
        %   0= don't print overall volumetric budget
        %   <>0 print overall volumetric budget at given frequency
        % ICBCFL
        %   0= don't write cell to cell flows to any file
        %  <>0 save cell to cell flows at given frequency to budget file depending
        %  on flag in the other package parameters I....CB

        % printing desired only if ihddfl>0
        % we use k as a print-frequency counter (within stress
        % period) and always save at the end of each stress period.

        if oc.INCODE(iP)==0  % Prevent>0 to treat all layers equally
            %3 Hdpr Ddpr Hdsv Ddsv    This is in fact layer specific
            fprintf(fid,'%10d%10d%10d%10d     Hdpr(%d) Ddpr(%d) Hdsv(%d) Ddsv(%d)\n',...
                oc.Hdpr(iP),oc.Ddpr(iP),oc.Hdsv(iP),oc.Ddsv(iP),k,k,k,k);
            %  Hdpr  <>0 prints head otherwise no printing
            %  Ddpr  <>0 prints ddn  otherwise no printing
            %  Hdsv  <>0 saves  head otherwise no saving
            %  Ddsv  <>0 saves  ddn  otherwise no saving
        end
    end
end

fclose(fid);
