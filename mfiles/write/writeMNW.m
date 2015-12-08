function writeMNW(basename,mnw)
%WRITEMNW writes input file for MODFLOW's NMW2 package
%
% writeMNW(basename,mnw)
%
% TO 110807

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if isempty(mnw.MNW)  % check if any data are provided for this package
    fprintf(['No data are provided for package <<%s>> !\n',...
        'Provide data or switch off this package in the NAM worksheet.\n'],...
        mnw.type);
end

fid=fopen([basename, '.',mnw.ext],'wt');

AsPrevious=-1; % flag indicating to reuse BCN of previous stress period

%% 0
fprintf(fid,'# MATLAB writeMNW for %s,  %s\n','MNW2 package',datestr(now));
fprintf(    '# MATLAB writeMNW for %s,  %s\n','MNW2 package',datestr(now));

%% 2

fprintf(fid,'%14d%14d             ; MAXACTMNW  MNWPRNT\n',length(mnw.MNW),mnw.MNWPRNT);
 
%% The non parameter (ordinary) values for this boundary condition

info='WELLID  Q';

mnwperlist = MNWprint(fid,mnw.MNW); % well data printed, stress data returned

for iPer=1:mnw.NPER   %length(mnw.ITMP)
   I2=find(mnwperlist(:,1)==iPer);  % which BCN on in this period?
   
   if iPer>1 && all(all(mnwperlist(I2,2:end)==mnwperlist(I1,2:end)))     
                fprintf(fid,'%14d                           ; AsPrevious       iPer(%d)\n',AsPrevious,iPer);
   else
       fprintf(fid,'%14d                           ; ITMP(%d). Following lines: %s\n',length(I2),iPer,info);
       for iw=1:length(I2)
           fprintf(fid,'%-20s %12g        ; WELLID(%d)   Q    iPER(%d)\n',...
               mnw.MNW(mnwperlist(I2(iw),2)).WELLID,...
               mnwperlist(I2(iw),3),...
               iw,iPer);
       end
   end
   I1=I2;
end

fclose(fid);
end

function mnwperlist=MNWprint(fid,MNW)
% MNWprint writes MNW well data

    %  info='Qdes MN QWval RW Skin Hlim Href DD LWGRP';
    %% NMW2 requires printing MNW info first
    for imnw=1:length(MNW)
        mnw=MNW(imnw);
        fprintf(fid,'%-20s %7d             ; WELLID(%d) NNODES\n',mnw.WELLID,mnw.NNODES,imnw);
        fprintf(fid,'%-24s %3d %3d %3d %3d ; LOSSTYPE PUMPLOC QLIMIT PPFLAG PUMPCAP\n',...
            mnw.LOSSTYPE,mnw.PUMPLOC,mnw.QLIMIT,mnw.PPFLAG,mnw.PUMPCAP);
        fprintf(fid,'  %12g %12g %12g ; rw RSKIN KSKIN\n',mnw.rw,mnw.rskin,mnw.kskin);
        for innodes=1:abs(mnw.NNODES)
            fprintf(fid,'  %12g %12g  %5d %5d ; zscrtop(%d) zscrbot(%d) ROW COL\n',...
                mnw.zscrtop(innodes),mnw.zscrbot(innodes),mnw.iy,mnw.ix,innodes,innodes);
        end
    end
    
    NPER=MNW(1).NPER;
    % The transient data
    mnwperlist=NaN(mnw.NPER*length(MNW),3);
    for iw=1:length(MNW)
        mnwperlist((iw-1)*NPER+(1:NPER),:) =...
            [(1:NPER)',ones(NPER,1)*iw,MNW(iw).Q]; % [stresperiod wellnr Q]
    end
    mnwperlist=sortrows(reshape(mnwperlist,[NPER*length(MNW),3]));
end
