function ZTA = cutoutZTA(ZTA,Ix,Iy)
%CUTOUTZTA cuts out a piece of ZTA accoring to input rows and columns
%
% Example:
%    ZTA = cutoutZTA(ZTA,Ix,Iy)
%
% Cutout ZTA(output of SWI package of same structure as BUD when read in by readBud([basename '.ZTA'])
% to reflect the new coordinate indices Ix and Iy
%
% These function are superseded by methods of the gridObj
%
% see also: gridObj modelObj cutoutBCN cutoutXTA RefineBCN RefindGrid removeCBD
%
% TO 120507

for it=1:length(ZTA)
    for iterm=1:length(ZTA(it).term)
        ZTA(it).term{iterm}=ZTA(it).term{iterm}(Iy,Ix,:);
        ZTA(it).NROW = length(Iy);
        ZTA(it).NCOL = length(Ix);
        ZTA(it).rows = 1:length(Iy);
        ZTA(it).cols = 1:length(Ix);
    end
end

