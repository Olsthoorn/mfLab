function writeLAK(basename,lak,transient)
%WRITELAK writes input file for MODFLOW's lake package
%
% Example:
%    writeLAK(basename,lak) -- write modflow lak file
%
% TO 070702


% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fid=fopen([basename, '.lak'],'wt');

%0
fprintf(fid,'# MATLAB writeLAK %s\r\n',datestr(now));
fprintf(    '# MATLAB writeLAK %s\n'  ,datestr(now));

%1
fprintf(fid,'%10i%10i\r\n',lak.NLAKES,lak.ILKCB);
%2
fprintf(fid,'%10.4f%10i%10.4f\r\n',lak.THETA,lak.NSSITR,lak.SSNCR);
%3
fprintf(fid,'%10.4f%10.4f%10.4f\r\n',lak.STAGES,lak.SSMN(1),lak.SSMX(1));
%3
for i=1:length(lak.ITMPLK)
    %4
    fprintf(fid,'%10i%10i%10i\r\n',lak.ITMPLK(i),lak.ITMP1(i),lak.LWRT(i));
    if lak.ITMPLK(i)>0
        %5
        fprintf(fid,'%10f\r\n',lak.LKARR(i));  % under construction
        %6
        fprintf(fid,'%10f\r\n',lak.BDLKNC(i));  % under construction
        %7
        fprintf(fid,'%5i\r\n',lak.NSLMS(i));
    end
    if lak.ITMPLK(i)>0 && lak.NSLMS(i)>0
        %8a and 8b
        fprintf(fid,'Use lak.ITMPLK and lak.NSLMS without subsystems illegal\r\n');
    end
    if lak.ITMP1(i)>=0
        %9a
        fprintf(fid,'%10.4f%10.4f%10.4f%10.4f',lak.PRCPLK(i),lak.EVAPLK(i),lak.RNF(i),lak.WTHDRW(i));
        if ~transient
            fprintf(fid,'%10.4f%10.4f',lak.SSMN(i),lak.SSMX(i));
        end
        fprintf(fid,'\r\n');
    end
    %9b and 9c not needed
end

fclose(fid);
