function  writeADV(basename,adv)
%WRITEADV writes basic advection transport package file
%
% Example:
%    writeADV(basename,adv)
%
% TO 0706030 081227

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

adv.unit=100;  % incompatibility with UD2REL (MT3DMS p97)

fid=fopen([basename,'.',adv.ext],'wt');

%HEADING 1+2 (<=80 chars)  -- must not be included in ADV file
fprintf(    '%s\n',['# MT3DMS writeADV ' datestr(now)]);


%B1 MIXELM PERCEL MXPART NADVFD (I10 F10.0 I10 I100
%   MIXELM adv solution option flag
%   PERCEL courant number, use 1
%   MXPART max total number of moving particles
%   NADVFD weighting scheme flag
fprintf(fid,'%10d%10g%10d%10d    MIXELM PERCEL MXPART NADVFD\n',...
    adv.MIXELM,adv.PERCEL,adv.MXPART,adv.NADVFD);

%B2 if MIXELM=1,2 or 3: ITRACK WD (I10 F10.0)
%   ITRACK particle tracking algorithm flag (p115)
%   WD     concentration weighting factor   (0.5 generally adequate)
if adv.MIXELM==1 || adv.MIXELM==2 || adv.MIXELM==3
    fprintf(fid,'%10d%10g     ITRACK WD\n',adv.ITRACK,adv.WD);
end

%B3 if MIXELM=1 or 3
%   DCEPS NPLANE NPL NPH NPMIN NPMAX (F10.0 5I10)
%   DCEPS  small relatie cell conce gradient
%   NPLANE random of fixed particle placement pattern flag 0=random
%   NPL initial # of particles in DCEPS cells, use 0
%   NPH initial # of particles in all other cells (16 in 2D to 32 in 3D adequate)
%   NPMIN min # of particles per cell use 2 (see p116) 
%   NPMAX max # of particles per cell use 2*NPH (see p116)
if adv.MIXELM==1 || adv.MIXELM==3
    fprintf(fid,'%10g%10d%10d%10d%10d%10d     DCEPS NPLANE NPL NPH NPMIN NPMAX\n',...
        adv.DCEPS,adv.NPLANE,adv.NPL,adv.NPH,adv.NPMIN,adv.NPMAX);
end

%B4 (if MIXELM=2 or 3) INTERP NLSINK NPSINK (3I10)
if adv.MIXELM==2 || adv.MIXELM==3
    adv.INTERP=1;            %   INTERP conc interpolation method flag, currently must be 1
    adv.NLSINK=adv.NPLANE;   %   NLSINK random or fixec particle placement flag in MMOC scheme use NLPLANE
    adv.NPSINK=adv.NPH;      %   NPSINK # particles to approximate sink cells in MMOC, use NPH
    fprintf(fid,'%10d%10d%10d     INTERP NLSINK NPSINK\n',...
        adv.INTERP,adv.NLSINK,adv.NPSINK);
end

%B5 (if MIXELM=3) DCHMOC (F10.0)
% DCHMOC crit rel conc grad controlling selective use of MOC or MMOC in HMOC
if adv.MIXELM==3
    fprintf(fid,'%10g      DCHMOC\n',adv.DCHMOC);
end

fclose(fid);
