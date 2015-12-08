function hfb=readHFB(fname,hfb)
%READHFB reads MODFLOW's horizontal flow barrier package
%
% Example:
% readHFB(basename,hfb);
%
% TO 070630 090713 120322

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fprintf('# MATLAB readHFB %s\n',datestr(now));

fid=fopen(fname,'r');

skipmodflowcomments(fid)

%1.
s=fgets(fid); %  NPHFB    MXFB     NHFBNP
C=textscan(s,'%f %f %f',1);
hfb.NPHFB  =C{1};
hfb.MXFB   =C{2};
hfb.NHFBNP =C{3};

%2  PARNAM   PARTYP   Parval   NLST
if hfb.NPHFB>0
    for iPar=1:hfb.NPHFB
        C=textscan(s,'%s %s %f %d',1); %    PARNAM   PARTYP   Parval   NLST
        hfb.PARNAM{iPar}=C{1};
        hfb.PARTYP{iPar}=C{2};
        hfb.Parval(iPar)=C{3};
        hfb.NLST(iPar)  =C{4};
        hfb.LST(iPar)=textscan(fid,'%f %f %f %f %f %f',hfb.NLST(iPar),'CollectOutput',1);
        %    Layer    IROW1    ICOL1    IROW2   ICOL2   Factor
    end
end

for iper=1:hfb.NPER
    %4
    if hfb.NHFBNP>0 
        hfb.cel(iper,1)=textscan(fid,'%f %f %f %f %f %f',hfb.NHFBNP,'CollectOutput',1);
        % Layer IROW1 ICOL1 IROW2 ICOL2 Hydchr
    end
    for iPar=1:hfb.NPHFB
        %5
        hfb.NACTHFB=fscanf(fid,'%d',[1,1]); fgets(fid);
        if hfb.NACTHFB>0
            %6
            for iN=1:hfb.NACTHFB
                hfb.Pname{iPar}=fgets(fid);
            end
        end
    end
end

fclose(fid);
