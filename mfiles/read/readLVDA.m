function  lvda=readLVDA(fname,lvda)
%READLVDA reads Layer Variabl Direction Anisotropy package file
%
% Example:
%    lvda=readLVDA(basename,lvda)
%
% TO 0706030 090713 090717

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

%0.
fprintf('# MATLAB readLVDA %s\n',datestr(now));

fid=fopen(fname,'r');
skipmodflowcomments(fid);

%1.
lvda.NPLVDA =fscanf(fid,'%d',1);  % number of lvda parameters
fprintf(fgets(fid));

for iPar=1:lvda.NPLVDA
    lvda.PARNAM{iPar}=fscanf(fid,'%s',1);
    lvda.PARTYP{iPar}=fscanf(fid,'%s',1);
    lvda.Parval(iPar)=fscanf(fid,'%f',1);
    lvda.NCLU(iPar)  =fscanf(fid,'%d',1);
    fprintf(fgets(fid));
%9  parameter clusters
    for iClu=1:lvda.NCLU(iPar)
        lvda.Cluster(iPar).Layer(iClu) =fscanf(fid,'%d',1);
        lvda.Cluster(iPar).Mltarr{iClu}=fscanf(fid,'%s',1);
        lvda.Cluster(iPar).Zonarr{iClu}=fscanf(fid,'%s',1);
        s=fgets(fid);
        lvda.Cluster(iPar).IZ{iClu}    =sscanf(s,'%d');
    end
end

fclose(fid);
