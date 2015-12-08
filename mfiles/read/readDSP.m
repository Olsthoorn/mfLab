function  readDSP(basename)
%READDSP reads MT3DMS's dispersion diffusion package input file
%
% Example:
%    readDSP(basename,dsp);
%
% TO 110112

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

dsp.unit=100;  % incompatibility with UD2REL (MT3DMS p97)

fid=fopen([basename,'.',dsp.ext],'wt');

%C1 HEADING No header allowed in DSP file
fprintf(    '%s\n',['# MT3DMS readDSP ' datestr(now)]);

%C0
if dsp.MultiDiffusion
    fprintf(fid,'$ MultiDiffusion\n');
end

%C1 AL, longitudinal dispersivity this is per CELL LAYERWISE
for i=1:dsp.NLAY
   dsp.AL=rarray(fid,[NROW,NCOL]); % long disp);
end

%C2 TRPT, aL/aTH, use >=0.1  THIS is 1 value per LAYER
dsp.TRPT=rarray(fid,[NLAY,1],'noheader');

%C3 TRPV, aL/aTV, use >=0.01 THIS is 1 value per LAYER
dsp.TRPV=rarray(fid,[NLAY,1],'noheader');

if dsp.MultiDiffusion
    for iComp=1:dsp.NCOMP
        dsp.DMCOEF{iComp}=NaN(NROW,NCOL,NLAY);
        for iLay=1:dsp.NLAY
            dsp.DMCOEF{iComp}(:,:,iLay)=rarray(fid,[NROW,NCOL]);
        end
    end
else
    for iComp=1:dsp.NCOMP
        %C4 DMCOEF,  >=0 or use 0 if insignificant (one value per layer)
        dsp.DMCOEF{iComp}=rarray(fid,[NLAY,1]);
    end
end

fclose(fid);
