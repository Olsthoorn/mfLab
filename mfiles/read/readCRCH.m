function crch=readCRCH(fname,pth,crch)
%READCRCH reads conduit flow recharge package input file
%
% Example:
%    crch=readCRCH([pth fname] ,cfp);
%
% TO 090708 090713 090714

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

%0
fprintf('# MATLAB readCRCH %s\n',datestr(now));
fid=fopen([pth fname],'r');

crch.IFLAG_CRCH =zeros(crch.NPER,1);
crch.P_CRCH=[];

for iPer=1:crch.NPER
    %1
    fprintf(fgets(fid));  % reaquired commentline
    crch.IFLAG_CRCH(iPer)=fscanf(fid,'%d',1); fgets(fid);

    %2
    if crch.IFLAG_CRCH(iPer)~=-1
        crch.P_CRCH(iPer).CRCH=zeros(crch.NNODES,2);
        for i=1:crch.NNODES
            crch.P_CRCH(iPer).CRCH(i,:)=fscanf(fid,'%f %f',[2,1]); fgets(fid);
        end
    end
end
fclose(fid);
