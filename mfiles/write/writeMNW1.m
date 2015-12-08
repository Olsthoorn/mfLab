function MAXNMW = writeMNW1(basename,mnw)
%WRITEMNW1 writs input file for multinode well (version 1) package
%
% Example:
%    writeMNW1(basename,mnw) -- write modflow MNW1 package file
%
% TO 110807

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if isempty(mnw.MNW)  % check if any data are provided for this package
    fprintf(['No data are provided for package <<%s>> !\n',...
        'Provide data or switch off this package in the NAM worksheet.\n'],...
        mnw.type);
end

fid=fopen([basename,'.',mnw.ext],'wt');

%% 0
fprintf(fid,'# MATLAB writeMNW1 for %s,  %s\n','MNW1 package',datestr(now));
fprintf(    '# MATLAB writeMNW1 for %s,  %s\n','MNW1 package',datestr(now));

%% total number of MNW cells

MAXNMW=0; for i=1:numel(mnw.MNW), MAXNMW=MAXNMW + numel(mnw.MNW(i).idx); end

mnw.IWELPT=0; % always print output

%% Writing records per simulation
%% #1
fprintf(fid,'%10d%10d%10d             !! MAXNMW ICB IWELPT\n',MAXNMW,mnw.ICB,mnw.IWELPT);

%% #2  Must be the same for all wells so take values of first well
fprintf(fid,'%-10s%10g                       !! LOSSTYPE (PLossMNW)\n',mnw.LOSSTYPE,mnw.MNW(1).lossP);
fprintf('LOSSTYPE will be %d and lossP = %g, taken from first MNW\n',mnw.LOSSTYPE,mnw.MNW(1).lossP);

%% #3
fu=freeunits(3,mnw.allUNITS);

% #3a
if mnw.WELL1
    fprintf(fid,'FILE:%-15s   WEL1:%-4d\n'           ,[basename '.MNWa'],fu(1));
end

% #3b
if mnw.BYNODE
    fprintf(fid,'FILE:%-15s   BYNODE:%-4d   ALLTIME\n',[basename '.MNWb'],fu(2));
end

% #3c
if mnw.QSUM
    fprintf(fid,'FILE:%-15s   QSUM:%-4d   ALLTIME\n',[basename '.MNWc'],fu(3));
end

%% The non parameter (ordinary) values for this boundary condition

% Asprevious doesn't work here, would make this routine unnaturally
% complicated
   
   %% #5
    for iper=1:mnw.NPER

       fprintf(fid,'# ===================================   S T R E S S   P E R I O D     =   %3d   ============================\n',iper);

       %% #4
       fprintf(fid,'%10d       !! ITMP\n',MAXNMW);
       fprintf(fid,...
           '#        L_________R_________C______Qdes___MN_DD_____QWVal________RW______SKIN______Hlim______Href____Iwgrp\n');

       mnw.MNW.write(fid,iper);
    end

fclose(fid);

end

function u=freeunits(n,UNIT)
    I=zeros(1,max(UNIT)+n); % total range of units used
    I(UNIT)=1;              % in use
    I(1:10)=1;              % assume also in use
    u=find(I==0);           % unused
    u=u(1:n);               % give me the first n of them
end

