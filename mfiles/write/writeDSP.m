function  writeDSP(basename,dsp)
%WRITEDSP writes input file for MODFLOW's dispersion diffusion package
%
% Example:
%    writeDSP(basename,dsp)
%
% TO 0706030 081227

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

dsp.unit=100;  % incompatibility with UD2REL (MT3DMS p97)

fid=fopen([basename,'.',dsp.ext],'wt');

%C0
if dsp.MultiDiffusion
    fprintf(fid,'$ MultiDiffusion\n');
end

%C1 HEADING 1+2 (<=80 chars) --- No header allowed in DSP file
fprintf(    '%s\n',['# MT3DMS writeDSP ' datestr(now)]);


%C1 AL, longitudinal dispersivity this is per CELL LAYERWISE
for i=1:dsp.NLAY
   if iscell(dsp.AL)
       warray(fid,dsp.AL{i},dsp.unit,'(10E12.4)',sprintf('AL(%d) (long disp)',i));
   else
       warray(fid,dsp.AL(:,:,i),dsp.unit,'(10E12.4)',sprintf('AL(%d) (long disp)',i));
   end
end

%C2 TRPT, aL/aTH, use >=0.1  THIS is 1 value per LAYER
warray(fid,dsp.TRPT(:)',dsp.unit,'(10E12.4)',sprintf('TRPT(1..%d) aTH/aL',dsp.NLAY));

%C3 TRPV, aL/aTV, use >=0.01 THIS is 1 value per LAYER
warray(fid,dsp.TRPV(:)',dsp.unit,'(10E12.4)',sprintf('TRPV(1..%d) aTV/aL',dsp.NLAY));

if dsp.MultiDiffusion
    if isa(dsp.DMCOEF,'cell')
        for iComp=1:dsp.NCOMP
            for iLay=1:dsp.NLAY
                warray(fid,dsp.DMCOEF{iComp}(:,:,iLay),dsp.unit,'(10E12.4)',sprintf('DMCOEF(iLay=%d, Species=%d)',iLay,iComp));
            end
        end
    else
       error(['%s: DMCOEF must be a cell array with one cell per species,\n',...
           'which holds a full 3D array with the cell values for that species.'],mfilename);
%        
%        % it's an array of NROW,NCOL,NLAY,NCOMP which is sometimes easier to specify in ma_adapt
%         for iComp=1:dsp.NCOMP
%             for iLay=1:dsp.NLAY
%                 warray(fid,squeeze(dsp.DMCOEF(:,:,iLay,iComp)),dsp.unit,'(10E12.4)',sprintf('DMCOEF(iLay=%d, Species=%d)',iLay,iComp));
%             end
%         end
    end
else
    %C4 DMCOEF,  >=0 or use 0 if insignificant (one value per layer)
    if isa(dsp.DMCOEF,'cell')
        for iComp=1:dsp.NCOMP
            for iLay = 1:dsp.NLAY
                warray(fid,dsp.DMCOEF{iComp}(iLay),dsp.unit,'(10E12.4)',sprintf('DMCOEF(Lay=1..%d) Species %d',dsp.NLAY,iComp));
            end
        end
    else
        error(['%s: DMCOEF must be a cell array with one cell per species containing a column that\n',...
        'holds one value per layer'],mfilename);
%         for iComp=1:dsp.NCOMP
%             for iLay = 1:dsp.NLAY
%                 warray(fid,dsp.DMCOEF{iComp}(iLay),dsp.unit,'(10E12.4)',sprintf('DMCOEF(Lay=1..%d) Species %d',dsp.NLAY,iComp));
%             end
%         end
    end
end
fclose(fid);
