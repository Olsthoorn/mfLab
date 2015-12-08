function  writeRCT(basename,rct)
%WRITERCT writes input file for MT3DMS's reaction package (RCT)
%
% Example:
%    writeRCT(basename,rct) --- write RCT chemical reaction package file
%
% TO 0706030 081227, 091109

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

SP1txt={'',         'Kd',         'Kf',         'Kl',         'Kd',      'dummy',         'Kd'};
SP2txt={'',      'dummy',        'a_f','C_SorbSites','massTrfRate','massTrfRate','massTrfRate'};
RC1txt={'','rr_dissSpec','rr_dissSpec','rr_dissSpec','rr_dissSpec','rr_dissSpec','rr_dissSpec'};
RC2txt={'','rr_sorbSpec','rr_sorbSpec','rr_sorbSpec','rr_sorbSpec','rr_sorbSpec','rr_sorbSpec'};

if rct.ISOTHM==0 && rct.IREACT~=0
    error(['\nError: You combined ISOTHM==0 with IREACT<>0.\n',...
        'This causes MT3DMS to produce wrong results.\n',...
        'This is a surely a bug in the MT3DMS code.\n',...
        'In your favor I will now stop the program!\n',...
        'As a workaround set ISOTHM=1 and SP1==0 and SP2==0.\n',...
        '\n',...
        'Sincerely yours\n',...
        'Theo Olsthoorn %s\n'],datestr(now));
end

fid=fopen([basename,'.RCT'],'wt');

%E1 HEADING 1+2 (<=80 chars)
%fprintf(fid,'%s\n',['# MATLAB writeRCT ' datestr(now)]); permitted in SWT but not in MT3DMS
 fprintf(    '%s\n',['# MT3DMS writeRCT ' datestr(now)]);

%E1 ISOTHM IREACT IRCTOP IGETSC (4I10)
%   ISOTHM adsorpton type flag
%   0=no sorption
%   1=linear sorption
%   2=freundlich
%   3=lamgmuir
%   4=first-order kinetic sorption
%   5=dual domain mass transfer (without sorption)
%   6=dual domain mass transfer (with sorption)
if rct.ISOTHM<0 || rct.ISOTHM>6,
    errstr=['ISOTHM in rct package equals %d, but must be between 1 and 6,...\n',...
           'Switch off RCT package in NAM sheet if no reactions or sorption is required !'];
    error(errstr,rct.ISOTHM);
end 

% IREACT=0 non kinetic reaction 1 first order irreversible reaction
% IRCTOP raction variables input method flag
%    >=2 3D reaction variable are input as 3D arrays
%    < 2 input as 1D array with each value in the arry correspondin to a
%    single layers (old)
% IGETSC flag to indicate input of initial sorbed concentrations these are
%    only read of IGETSC=1 in case ISOTHM=4, 5 or 6
fprintf(fid,'%10d%10d%10d%10d     ISOTHM IREACT IRCTOP IGETSC\n',...
    rct.ISOTHM,rct.IREACT,rct.IRCTOP,rct.IGETSC);

%E2A Enter if ISOTHM=1, 2, 3, 4 or 6 but not 5 RHOB(NROW,NCOL);
if rct.ISOTHM~=5
    for iLay=1:rct.NLAY
        warray(fid,rct.RHOB(:,:,iLay),rct.unit,'(10E15.6)',sprintf('RHOB{%d}',iLay),true,rct.FREE);
    end
end

%E2B, if ISOTHM=5 or 6, PRSITY2(NCOL,NROW)
if rct.ISOTHM==5 || rct.ISOTHM==6
   for iLay=1:rct.NLAY
       warray(fid,rct.PRSITY2(:,:,iLay),rct.unit,'(10E15.6)',sprintf('PRSITY2{%d}',iLay),true,rct.FREE);
   end
end

%E2C if IGETSC>0 initial sorbed conc for each species
if rct.ISOTHM>0 && rct.IGETSC>0
    for iComp=1:rct.NCOMP
        for iLay=1:rct.NLAY
           warray(fid,rct.SRCONC{iCOMP}(:,:,iLay),rct.unit,'(10E15.6)',...
               sprintf('SRCONC{%d} layer(%d)',iComp,iLay),true,rct.FREE);
        end
    end
end

%E3 for each species if ISOTHM>0 (sorption always)
if rct.ISOTHM>0
    for iComp=1:rct.NCOMP
        for iLay=1:rct.NLAY
            warray(fid,rct.SP1{iComp}(:,:,iLay),rct.unit,'(10E15.6)',...
                sprintf('SP1: %s{%d} layer(%d)',SP1txt{rct.ISOTHM+1},iComp,iLay),true,rct.FREE);
        end
    end
end

%E4 for each species if ISOTHM>0 (sorption active)
if rct.ISOTHM>0
    for iComp=1:rct.NCOMP
        for iLay=1:rct.NLAY
            warray(fid,rct.SP2{iComp}(:,:,iLay),rct.unit,'(10E15.6)',...
                sprintf('SP2: %s{%d} layer(%d)',SP2txt{rct.ISOTHM+1},iComp,iLay),true,rct.FREE);
        end
    end
end

%E5 for each species if REACT>0 (decay active)
if rct.IREACT>0
    for iComp=1:rct.NCOMP
        for iLay=1:rct.NLAY
            warray(fid,rct.RC1{iComp}(:,:,iLay),rct.unit,...
                '(10E15.6)',sprintf('RC1: %s{%d} layer(%d)',RC1txt{rct.ISOTHM+1},iComp,iLay),true,rct.FREE);
        end
    end
end

%E6 for each species if IREACT>0 (decay active)
if rct.IREACT>0
    for iComp=1:rct.NCOMP
        for iLay=1:rct.NLAY
            warray(fid,rct.RC2{iComp}(:,:,iLay),rct.unit,...
                '(10E15.6)',sprintf('RC2: %s{%d} layer(%d)',RC2txt{rct.ISOTHM+1},iComp,iLay),true,rct.FREE);
        end
    end
end

fclose(fid);
