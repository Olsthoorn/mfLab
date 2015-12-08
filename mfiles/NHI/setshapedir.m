function setshapedir(pth,basename)
%SETSHAPEDIR renames files in shapefile dir so that they all have the same basenae and
%
% Example:
%    setshapedir(pth,basename)  (no output)
%
% proper extension
%
% TO 120401

olddir=pwd;
cd(pth)
d=dir([basename '*']);
for iFile=1:length(d)
    if d(iFile).name(1)~='.'
        i=findstr('.',d(iFile).name);
        if isempty(i)
            switch d(iFile).name(end-2:end)
                case {'dbf','prj','qix','sbn','sbx','shp','shx'}
                    fprintf('Renaming ''%s'' to ''%s''\n',[pth '\' d(iFile).name],[d(iFile).name(1:end-3) '.' d(iFile).name(end-2:end)]);
                    [SUCCESS,MESSAGE]=movefile([pth '\' d(iFile).name],[d(iFile).name(1:end-3) '.' d(iFile).name(end-2:end)]);
                    if ~SUCCESS,
                        error(MESSAGE);
                    end
            end
        end
    end
end

cd(olddir);