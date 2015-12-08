function writeOCwords(basename,oc)
%WRITEOCWORDS writes output control input file using WORDS mode
%
% Example:
%    writeOC(fname,oc)
%
% See manual of MODFLOW 2000, the Open File Report 00-92, p52-53
%
% TO 130204
%
% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fid=fopen([basename, '.',oc.ext],'wt');

%0
fprintf(fid,'# MATLAB %s %s\n',mfilename,datestr(now));
fprintf(    '# MATLAB %s %s\n',mfilename,datestr(now));

%1 IHEDFM IDDNFM IHEDUN IDDNUN (For each simulation)
fprintf(fid,'HEAD PRINT FORMAT %d\n',    oc.IHEDFM);
fprintf(fid,'HEAD SAVE UNIT %d\n',       oc.IHEDUN);
fprintf(fid,'DRAWDOWN PRINT FORMAT %d\n',oc.IDDNFM);
fprintf(fid,'DRAWDOWN SAVE UNIT %d\n',   oc.IDDNUN);
if oc.COMPACT
    fprintf(fid,'COMPACT BUDGET FILES\n');
end

%2
for iP=1:oc.NPER  % for each stress period

    k=0; % in-stress-period time-step counter; reset at beginning of each stress period
    
    for iS=1:oc.NSTP(iP)  % for all time steps within this stress period

        k=k+1; % counting of total number of time steps
        
        ihddfl=oc.IHDDFL(iP) & (rem(k,oc.IHDDFL(iP))==0 | k==oc.NSTP(iP));
        ibudfl=oc.IBUDFL(iP) & (rem(k,oc.IBUDFL(iP))==0 | k==oc.NSTP(iP));
        icbcfl=oc.ICBCFL(iP) & (rem(k,oc.ICBCFL(iP))==0 | k==oc.NSTP(iP));
        
        %2 INCODE IHDDFL IBUDSV ICBCFL
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

        if ihddfl || ibudfl || icbcfl
            fprintf(fid,'\nPERIOD %d STEP %d\n',iP,iS);
            
            if ihddfl
                if oc.Hdpr(iP), fprintf(fid,'PRINT HEAD\n'); end
                if oc.Hdsv(iP), fprintf(fid,'SAVE HEAD\n');  end
                if oc.Ddpr(iP), fprintf(fid,'PRINT DRAWDOWN\n'); end
                if oc.Ddsv(iP), fprintf(fid,'SAVE DRAWDOWN\n');  end
            end
            
            if ibudfl
                fprintf(fid,'PRINT BUDGET\n');
            end
            if icbcfl
                fprintf(fid,'SAVE BUDGET\n');
            end
        end
    end
end

fclose(fid);
