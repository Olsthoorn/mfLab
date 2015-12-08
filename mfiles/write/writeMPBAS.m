function writeMPBAS(basename,mpbas)
%WRITEMPBAS writes input file for MODPATH6's basic package
%
% Example:
%    writeMPBAS(basename,mpbas)
%
% TO 070630 130123

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fid=fopen([basename,'.',mpbas.ext],'wt');

recDesired = true;

%0.
fprintf(fid,'# MATLAB  writeMPBAS %s\n',datestr(now));
fprintf(    '# MODPATH writeMPBAS %s\n',datestr(now));

%1. value used for inactive cells
fprintf(fid,'%10g %10g     HNOFLO HDRY\n',mpbas.HNOFLO,mpbas.HDRY);

%2.
fprintf(fid,'%d  # of default iface labels\n',numel(find(~isnan(mpbas.defaultIFace))));
for iFace= find(~isnan(mpbas.defaultIFace))
    %3
    fprintf(fid,'%s\n',mpbas.budgetLabel{iFace});
    %4
    fprintf(fid,'%d\n',mpbas.defaultIFace(iFace));
end


%5. initial heads
warray(fid,mpbas.LAYTYP(:)'          ,mpbas.unit,'(40I2)','LAYTYP',~recDesired);

%6
for ilay=1:mpbas.NLAY
    warray(fid,mpbas.IBOUND(:,:,ilay),mpbas.unit,'(40I4)',sprintf('IBOUND{%d}',ilay),recDesired);
end

k=0;
for ilay=1:mpbas.NLAY
    %7
    warray(fid,mpbas.porosity(:,:,ilay),  mpbas.unit,'(10E12.3)',sprintf('porosity{%d}',ilay),recDesired);
    if mpbas.LAYCBD(ilay)
        k=k+1;
        %8
        warray(fid,mpbas.porosityCB(:,:,k),mpbas.unit,'(10E12.3)',sprintf('porosityCB{%d}',ilay),recDesired);
    end
end
fclose(fid);
